#    Glycogen_analysis.py (C) 2015, Heikki Lehvaeslaiho, Corrado Cali, Jumana Baghabra, Daniya Boges
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see http://www.gnu.org/licenses/

bl_info={
        "name":"Glycogen Analysis",
        "author":"Heikki Lehvaeslaiho, Corrado Cali, Jumana Baghabra, Daniya Boges",
        "version":"1.0",
        "location":"VIEW3D UI > Glycogen Analysis",
        "description":"Glycogen Analysis (clustering + Measurements + calculate Dbscan Optimum values)",
        "category":"objects"
}
import bpy
from bpy.props import*
import bmesh
import numpy as np
from numpy import linalg
from sklearn.cluster import DBSCAN
from sklearn import metrics
from random import random

import re
import csv
from scipy.stats import itemfreq
from collections import Counter
from collections import OrderedDict

from scipy.spatial import cKDTree
from scipy.spatial import distance
import time

import os
from os.path import expanduser
import time

# for OSX 10.8:
#export PYTHONPATH=/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages:$PYTHONPATH
# CentOS 7:
#export PYTHONPATH=/usr/local/lib/python3.4/site-packages:$PYTHONPATH && export PYTHONPATH=/usr/local/lib/python3.4:$PYTHONPATH && export PYTHONPATH=/usr/local/lib/python3.4/lib-dynload:$PYTHONPATH && export OMP_NUM_THREADS='20'
# SL6-vislab:
#export PYTHONPATH=/usr/local/lib/python3.4/lib-dynload:$PYTHONPATH	
#export OMP_NUM_THREADS='20'

#---------------------------------------------------------------------------------
# 	  					# TIMER #
#---------------------------------------------------------------------------------
class ModalTimerOperator(bpy.types.Operator):
    """Operator which runs its self from a timer"""
    bl_idname = "wm.modal_timer_operator"
    bl_label = "Modal Timer Operator"

    _timer = None

    def modal(self, context, event):
        if event.type in {'RIGHTMOUSE', 'ESC'}:
            return self.cancel(context)

        if event.type == 'TIMER':
            # change theme color, silly!
            color = context.user_preferences.themes[0].view_3d.space.gradients.high_gradient
            color.s = 1.0
            color.h += 0.01

        return {'PASS_THROUGH'}

    def execute(self, context):
        wm = context.window_manager
        self._timer = wm.event_timer_add(0.1, context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
#---------------------------------------------------------------------------------
# 	  					# INFO TEXTS #
#---------------------------------------------------------------------------------
class VIEW3D_OT_show_info(bpy.types.Operator):
	bl_idname = "view3d.show_info"
	bl_label = ""
	bl_description = "1- Autoset DBSCAN parameters | 2- Manual DBSCAN setup"
	
	def invoke(self,context,event):
		if not bpy.types.Scene.prop_dbscan_info:
			bpy.types.Scene.prop_dbscan_info = True
		else:
			bpy.types.Scene.prop_dbscan_info = False
			
		return{"FINISHED"}

#---------------------------------------------------------------------------------
# 	  OPERATOR ACTIVATE/RELOAD activates Selection, Clustering and Display Clusters
#---------------------------------------------------------------------------------
class VIEW3D_OT_activate_addon_button(bpy.types.Operator):
	bl_idname = "view3d.activate_addon_button"
	bl_label = " Activate/Reload addon"

	def func_init(self):
		bpy.types.Scene.data_names = []
		bpy.types.Scene.prop_obj_names = None
		bpy.types.Scene.data_obj_coords = None
		bpy.types.Scene.prop_bool_glycogen = None
		bpy.types.Scene.prop_bool_clusters = None

		bpy.types.Scene.prop_min_samples = None
		bpy.types.Scene.prop_eps1 = None
		bpy.types.Scene.prop_eps2 = None
		bpy.types.Scene.prop_interval = None
		#step2
		bpy.types.Scene.prop_min_samples_s2 = None #in sync with prop_min_samples step 1
		bpy.types.Scene.prop_optimum_eps = None
		bpy.types.Scene.prop_nclusters = None
		bpy.types.Scene.prop_silh = None
		#measurements: The Activate/Reload button if its pressed for Reload then it shouldnt get rid of certain existing values:
		if bpy.types.Scene.data_glyc_clusters is None:
			bpy.types.Scene.flag_clusters_measure = None
			bpy.types.Scene.flag_clusters_measured = None



	def func_load_objectsNames(self):
		lst = []
		symbols = "-_"
		for ob in bpy.data.objects:
			if ob.type != "MESH" and ob.type != "CURVE":
				continue
			ob.hide = False
			ob.select = False
			if str.isalpha(ob.name):
				lst.append(ob.name)
			else:
				for ch in ob.name:
					if str(ch).isdigit() or ch in symbols:
						common = ob.name[:ob.name.index(ch)]
						if common not in lst:
							lst.append(common)
						break
		if "Glycogen" in lst: #jan 10th/2015
			print("Created")
			bpy.types.Scene.prop_bool_glycogen = bpy.props.BoolProperty(name="Glycogens", description=" ",
			update=update_prop_bool_glyc)
			bpy.types.Scene.prop_bool_clusters = bpy.props.BoolProperty(name="Clusters", description=" ",
			update=update_prop_bool_clust)
			bpy.types.Scene.prop_min_samples = bpy.props.IntProperty(name='Minimum Samples', default=20)
			bpy.types.Scene.prop_eps1 = bpy.props.FloatProperty(name='from', default=0.2)
			bpy.types.Scene.prop_eps2 = bpy.props.FloatProperty(name='to', default=0.63)
			bpy.types.Scene.prop_interval = bpy.props.FloatProperty(name='Interval', default=0.01)
			#step 2
			bpy.types.Scene.prop_min_samples_s2 = bpy.props.IntProperty(name='Minimum Samples', default=20)
			bpy.types.Scene.prop_optimum_eps = bpy.props.FloatProperty(name='Optimum Epsilon', default=00)
			bpy.types.Scene.prop_nclusters = StringProperty(name='', default="")
			bpy.types.Scene.prop_silh = StringProperty(name='',default="")#Silhouette Coefficient
			if bpy.types.Scene.data_glyc_clusters is None:
				bpy.types.Scene.flag_clusters_measure = False #TEMP FOR TESTING GARFY
				bpy.types.Scene.flag_clusters_measured = False
			#glycogen Measures:
			empty = [("","","")]
			bpy.types.Scene.prop_glyc_neighbours = EnumProperty(name="Neibouring Objects",items=empty)
			bpy.types.Scene.prop_associated_glyco = EnumProperty(name='Associated Granules:',items=empty)
			bpy.types.Scene.prop_total_granules = StringProperty(name='Total Granules:',default="")
			bpy.types.Scene.prop_glyc_to_neighb_dist = StringProperty(name='Distance:',default="")
			bpy.types.Scene.data_glyc_distances = []
			if bpy.types.Scene.data_glyc_distances:
				print("its not none")

		else:
			if bpy.types.Scene.prop_bool_glycogen != None:
				del bpy.types.Scene.prop_bool_glycogen
			if bpy.types.Scene.prop_bool_clusters != None:
				del bpy.types.Scene.prop_bool_clusters
			#clustering option:
			if bpy.types.Scene.prop_min_samples:
				del bpy.types.prop_Scene.prop_min_samples
			if bpy.types.Scene.prop_eps1:
				del bpy.types.Scene.prop_eps1
			if bpy.types.Scene.prop_eps2:
				del bpy.types.Scene.prop_eps2
			if bpy.types.Scene.prop_interval:
				del bpy.types.Scene.prop_interval
			#glycogen option:
		return lst
	#---------------#
	def func_updateObjectsNames(self):
		lcl_list = self.func_load_objectsNames()
		objects = []
		bpy.types.Scene.prop_obj_names = []
		if lcl_list:
			for index, name in enumerate(lcl_list):
				objects.append((name, name, str(index))) #blender needs a tuple for a dropdown list UI
				bpy.context.scene.data_names.append(name)
			bpy.types.Scene.prop_obj_names = bpy.props.EnumProperty(items=objects)
	#---------------#
	def invoke(self,context,event):
		self.func_init()
		self.func_updateObjectsNames()
		return{"FINISHED"}
	#http://wiki.blender.org/index.php/Dev:2.5/Py/Scripts/Cookbook/Code_snippets/Interface#Invoke_versus_execute
	#def execute(self, context): #instaiate a new object with everychange/ new call
#---------------------------------------------------------------------------------
# 	  					# UPDATE / PUBLIC FUNCTIONS #
#---------------------------------------------------------------------------------
def update_prop_bool_glyc(self,context):
	#Toggle between two boolean values
	if bpy.context.scene.prop_bool_glycogen != bpy.context.scene.prop_bool_clusters:
		return None #to avoid infinite recursion
	
	if bpy.context.scene.prop_bool_glycogen == False:
		bpy.context.scene.prop_bool_clusters = True
		return None
	else:
		bpy.context.scene.prop_bool_clusters = False
		return None
	return None

def update_prop_bool_clust(self,context):
	if bpy.context.scene.prop_bool_clusters != bpy.context.scene.prop_bool_glycogen:	
		return None #continue

	if bpy.context.scene.prop_bool_clusters == False:
		bpy.context.scene.prop_bool_glycogen = True
		return None
	else:
		bpy.context.scene.prop_bool_glycogen = False
		return None
	return None

#called when neuro list selection is changed -it changes value of no.of.glycogens and associtaed glycogens list
def update_glyc_neighbours_dropdown(self,context):#updateFun1

	bpy.context.scene.prop_total_granules = str(bpy.types.Scene.data_noOfGlyc[int(bpy.context.scene.prop_glyc_neighbours)])
	
	listtemp = []
	listtemp = bpy.types.Scene.neuro_glyList_dict[bpy.types.Scene.data_glyc_neighbours[int(bpy.context.scene.prop_glyc_neighbours)]]
	#for associated gly, we need to update bpy.types.Scene not bpy.types.Context, as we want to change the entire dropdown list not an element of an existing one
	#updating associated granules
	enum2 = []
	for _index, enumval in enumerate(listtemp):
		enum2.append((enumval,enumval, ""))
	bpy.types.Scene.prop_associated_glyco = EnumProperty(name="Associated Granules:", items=enum2, update=update_associated_granules_dropdown)
	#updating distance
	for glyname, dist in bpy.context.scene.data_glyc_to_neighb_dist:
		if glyname == bpy.context.scene.prop_associated_glyco:
			bpy.context.scene.prop_glyc_to_neighb_dist = str(dist)

	return None

#associated granules update function = updates value of distances
def update_associated_granules_dropdown(self,context):
	for glyname, dist in bpy.context.scene.data_glyc_to_neighb_dist:
		if glyname == bpy.context.scene.prop_associated_glyco:
			bpy.context.scene.prop_glyc_to_neighb_dist = str(dist)

	return None

def updateFun3(self,context):

	return None

def updateFun4(self,context):

	return None

def func_unhide():
	for ob in bpy.data.objects:
		ob.hide = False
#---------------#
def func_unselect():
	for obj in bpy.data.objects:
		obj.select = False

def func_select(pattern):
	selected_objects=[]
	for ob in bpy.data.objects:
		if ob.type != 'MESH' and ob.type != 'CURVE':
			continue
		match = re.search(pattern+'*', ob.name)
		if not match:
			ob.hide = True
			continue
		ob.select = True
		selected_objects.append(ob)
	return selected_objects

def func_median_location(selected):
	#all objects are selected already, we just need to set workspace to the right layer
	for ob in selected:
		thisObjLayer = [i for i in range(len(ob.layers)) if ob.layers[i] == True] #returns objects layer
		if(thisObjLayer):
			layer_indx = thisObjLayer[0]
			bpy.context.scene.layers[layer_indx] = True #this will set workspace at THIS layer
			break
	bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')

def func_toGlobal_coords(ob): #NOT USED
		coord_matrix = ob.matrix_world
		for v in ob.data.vertices:
			loc = v.co
			globalCoords = coord_matrix * loc
			if ob.parent is None:
				bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
			else:
				bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
	
def populate_glycogen_coords():
	func_unhide()
	func_unselect()
	selected = func_select("Glycogen")
	bpy.types.Scene.data_obj_coords = []
	func_median_location(selected)

	for ob in selected:
		if ob.parent is None:
			bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(ob.location.x),str(ob.location.y),str(ob.location.z)])
		else:
			bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(ob.location.x),str(ob.location.y),str(ob.location.z)])

def glyc_toGlobal_coords(ob,objs_attrib,objs_verts):
    coord_matrix = ob.matrix_world
    for v in ob.data.vertices:
        loc = v.co
        globalCoords = coord_matrix *loc
        n = ob.name.rsplit('_',1)#why is it needed? str(n[0])
        if ob.parent is None:
            objs_attrib.append([n[0], "None"])
            objs_verts.append([str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
        else:
            objs_attrib.append([n[0],ob.parent.name])
            objs_verts.append([str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
    return(objs_attrib,objs_verts)

# getVertices will seperate vertices from attributes
def getVertices(pattern,coords_type): 
    func_unhide()
    func_unselect()
    selected_objects = []
    objs_attrib = [] #either name, parent or just name
    objs_verts = []
    objs_attrib_matrix = []
    objs_verts_matrix = []

    for ob in bpy.data.objects:
        if ob.type !='MESH':
            continue
        match = re.search(pattern+'*', ob.name)
        if not match:
            ob.hide = True
            continue
        else:
        	ob.select = True
        	selected_objects.append(ob)
 
    if selected_objects and coords_type == "Center Vertices":
        func_median_location(selected_objects)
        
        for ob in selected_objects:
        	#skipping the first granule as old script, we can leave that
            #if ob.name == 'Glycogen': 
            #    continue
            n = ob.name.rsplit('_',1)
            if ob.parent is None:
                objs_attrib.append([n[0],"None"])
                objs_verts.append([str(ob.location.x),str(ob.location.y),str(ob.location.z)])
            else:
                objs_attrib.append([n[0],ob.parent.name])
                objs_verts.append([str(ob.location.x),str(ob.location.y),str(ob.location.z)])

    elif selected_objects and coords_type == "All Vertices":
        for ob in selected_objects:
            objs_attrib,objs_verts =glyc_toGlobal_coords(ob,objs_attrib,objs_verts)
    
    return(np.array(objs_attrib),np.array(objs_verts))

def get_closest_distance(first_verts, second_verts):
	closests_points_info = []
	#first_verts = np.asarray(bpy.types.Scene.glycogen_verts_np)
	#second_verts = np.asarray(bpy.types.Scene.neur_obj_verts_np)
	first_verts = np.asarray(first_verts)
	second_verts = np.asarray(second_verts)
	
	#scipy cKDTree algorithm
	mytree = cKDTree(second_verts)
	dist, indexes = mytree.query(first_verts)
	dist= np.asarray(dist)
	indexes = np.asarray(indexes)
	a = np.vstack(dist)
	b = np.vstack(indexes)
	c = np.hstack((b,a))

	for i in range(0,len(c)):
		closests_points_info.append([i, c[i,0] , c[i,1]])
	
	return closests_points_info

#---------------------------------------------------------------------------------
# 	  								GLYCOGEN
#---------------------------------------------------------------------------------
class OBJECTS_OT_glycogens_nearest_neighbours(bpy.types.Operator):
	bl_idname= "objects.glycogens_neareste_neighbors"
	bl_label = ""
	
	def invoke(self,context,event):
		self.initialise()
		#bpy.ops.wm.modal_timer_operator()
		bpy.types.Scene.glycogen_attrib_np, bpy.types.Scene.glycogen_verts_np = getVertices("Glycogen", "Center Vertices")

		patterns = []
		if "bouton" in bpy.context.scene.data_names:
			patterns.append('bouton')
		elif "Bouton" in bpy.context.scene.data_names:
			patterns.append('Bouton')
		if "spine" in bpy.context.scene.data_names:
			patterns.append('spine')
		elif "Spine" in bpy.context.scene.data_names:
			patterns.append('Spine')

		firstTime = True

		for patt in patterns:
			temp_attrib_np, temp_verts_np = getVertices(patt, "All Vertices")
			if temp_attrib_np.size:
					print("getting vertices for: ", patt)
					if firstTime:
						bpy.types.Scene.neur_obj_verts_np = temp_verts_np
						bpy.types.Scene.neur_obj_attrib_np = temp_attrib_np
						firstTime = False
						#to set the big array to same dimensions of concatinated arrays
					else:
						bpy.types.Scene.neur_obj_attrib_np =np.concatenate((bpy.types.Scene.neur_obj_attrib_np, temp_attrib_np),axis=0)
						bpy.types.Scene.neur_obj_verts_np = np.concatenate((bpy.types.Scene.neur_obj_verts_np, temp_verts_np),axis=0)

		if bpy.types.Scene.neur_obj_attrib_np.size:
			bpy.types.Scene.data_glyc_distances = get_closest_distance(bpy.types.Scene.glycogen_verts_np, bpy.types.Scene.neur_obj_verts_np)
			bpy.types.Scene.data_glyc_distances_occur = self.occurences()


		# print
		#for i in range(0, len(bpy.types.Scene.data_glyc_distances)):
		#	print(bpy.types.Scene.data_glyc_distances[i])
			#if i == 5:
			#	break
		#for glyIndx, neuroIndx, dist in bpy.types.Scene.data_glyc_distances:
			#print(bpy.types.Scene.glycogen_attrib_np[glyIndx,0], bpy.types.Scene.neur_obj_attrib_np[neuroIndx,0], bpy.types.Scene.neur_obj_attrib_np[neuroIndx,1], dist)
		#print('NeuroObj','TotalAssociatedG','GlycogenList','distance', 'distance[0]')
		#for neuroObj, glyList in bpy.types.scene.neuro_glyList_dict.items():
		#	for glyName in glyList:
		#		gdistance = [dist for gname, dist in bpy.types.Scene.data_glyc_distances if gname == glyName]
		#		if gdistance:
		#			print(neurObj, len(glyList), glyName, gdistance, gdistance[0])
		self.visualization()
		return{"FINISHED"}

	def visualization(self):
		bpy.types.Scene.data_glyc_neighbours = []
		bpy.types.Scene.data_noOfGlyc = []
		bpy.types.Scene.data_glyc_to_neighb_dist = []
		enum1 = [] #neural strcure (dropdown)
		enum2 = [] #associated granules (dropdown)

		''' we need proper data strcuture for: neural strcure (dropdown), associated granules (dropdown),
 		NoGranules (Textbox), distance (Textbox)'''
		for neural_obj_name, noOfGlycogens in bpy.types.Scene.data_glyc_distances_occur:
			bpy.types.Scene.data_glyc_neighbours.append(neural_obj_name)
			bpy.types.Scene.data_noOfGlyc.append(noOfGlycogens)
			#1- neural neighbors (dropdown)
		for _index, enumval in enumerate(bpy.types.Scene.data_glyc_neighbours):
			enum1.append((str(_index), enumval, "")) ##the value in the middle will show in the UI
		bpy.types.Scene.prop_glyc_neighbours = EnumProperty(name="Neibouring Objects", items=enum1, 
			update=update_glyc_neighbours_dropdown)
		#2- associated granules (dropdown)
		# initilised while update - for now we will get the current value (1st row according to )
		for _index, enumval in enumerate(bpy.types.Scene.neuro_glyList_dict[
			bpy.types.Scene.data_glyc_neighbours[
			int(bpy.context.scene.prop_glyc_neighbours)]]):
			enum2.append((enumval,enumval, ""))
		#print("associated granules in visualizing_measurements:", enum2)
		bpy.types.Scene.prop_associated_glyco = EnumProperty(name='Associated Granules:', items=enum2, 
			update=update_associated_granules_dropdown)

		#3- Total Granules# (textbox)
		bpy.types.Scene.prop_total_granules = StringProperty(name='Total Granules:',
			default=str(bpy.types.Scene.data_noOfGlyc[int(bpy.context.scene.prop_glyc_neighbours)]),
			update=updateFun3)
		
		#4- Distance (textbox)
		for glyIndx, neuroIndx, dist in bpy.types.Scene.data_glyc_distances:
			bpy.types.Scene.data_glyc_to_neighb_dist.append(
				[bpy.types.Scene.glycogen_attrib_np[glyIndx,0],dist]
				)
		for glyname, dist in bpy.types.Scene.data_glyc_to_neighb_dist:
			if glyname == bpy.context.scene.prop_associated_glyco:
				bpy.types.Scene.prop_glyc_to_neighb_dist = StringProperty(name='Distance:',
					default=str(dist),
					update=updateFun4)

	def initialise(self):
		
		bpy.types.Scene.data_glyc_distances = []
		bpy.types.Scene.data_glyc_distances_occur = []

	def occurences(self):
		objects_names=[]
		#for the glycogen occurance output UIdropdownList:
		gly_names=[] 
		bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = {} #used for data export
		bpy.types.Scene.neuro_glyList_dict = {} #used for UI display (a group-by done on 'sorted')
		dict_temp = {}
		
		closest_points_np = (np.array(bpy.types.Scene.data_glyc_distances)).astype(int)
		
		for k in range(0, len(closest_points_np)):
			objects_names.append(" ".join
				((
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],0],
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],1]
				))
				)#join (object name, parent) with a " " between them
			gly_names.append(bpy.types.Scene.glycogen_attrib_np[closest_points_np[k,0],0])
		#end loop

		'''now we can create dictionary from obj&gly names:'''
		for i in range(0,len(gly_names)):
			#the below method in populating dictionaries will perform (itemfreq) on objects names so that it will not duplicate keys (robust?)
			dict_temp[gly_names[i]] = objects_names[i] 

		#sort the dictionary:b = OrderedDict(sorted(a.items()))
		bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = OrderedDict(sorted(dict_temp.items(), key=lambda val: val[1])) #checked
		
		#now we need to group glycogens in a list per neuro
		#===============================================
		#Populating bpy.types.Scene.neuro_glyList_dict
		templist1 = []
		current_neuro = list(bpy.types.Scene.neuro_gly_glyFreq_dic_sorted.values())[0]
		for glyname, objname in bpy.types.Scene.neuro_gly_glyFreq_dic_sorted.items(): #.item refers to a pair(key,value), switching keys to values and values to keys
			print("outer",current_neuro)
			if objname == str(current_neuro):
				print("in if",current_neuro)
				templist1.append(glyname)
			else:
				print("in else",current_neuro)
				bpy.types.Scene.neuro_glyList_dict[current_neuro]=templist1
				current_neuro = objname
				templist1 = []
				templist1.append(glyname)
		bpy.types.Scene.neuro_glyList_dict[current_neuro]=templist1 #store last list, as loop finishes before it gets to store it in'else'

		#===============================================

		#now we count instances of glycogens per neuro object closests to it
		templist=[]
		countInstances = Counter(objects_names) #countInstances is a dictionary
		#{'objName parent': itemfreq()...}
		for neur_obj_name, noOfGlycogens in countInstances.items():
			templist.append([neur_obj_name,noOfGlycogens])
		print("countInstances done")

		return templist
#--------------------------------------------------------------------------------
#							OPERATOR DISPLAY SELECTION
#--------------------------------------------------------------------------------
class VIEW3D_OT_display_selected(bpy.types.Operator):
	bl_idname = "view3d.display_selected"
	bl_label = "Display Selected"

	@classmethod
	def func_unhide(self):
		for ob in bpy.data.objects:
			ob.hide = False
	#---------------#
	def func_unselect(self):
		for obj in bpy.data.objects:
			obj.select = False
	#---------------#
	def func_toGlobal_coords(self,ob):
		coord_matrix = ob.matrix_world
		for v in ob.data.vertices:
			loc = v.co
			globalCoords = coord_matrix * loc
			if ob.parent is None:
				bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
			else:
				bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(globalCoords.x),str(globalCoords.y),str(globalCoords.z)])
	#---------------#
	def func_select(self,pattern):
		"""selects object with a match hit """
		selected_objects=[]
		for ob in bpy.data.objects:
			if ob.type != 'MESH' and ob.type != 'CURVE':
				continue
			match = re.search(pattern+'*', ob.name)
			if not match:
				ob.hide = True
				continue
			ob.select = True
			selected_objects.append(ob)
		return selected_objects
	#---------------#
	def func_median_location(self,selected):
		#all objects are selected already, we just need to set workspace to the right layer
		for ob in selected:
			thisObjLayer = [i for i in range(len(ob.layers)) if ob.layers[i] == True] #returns objects layer
			if(thisObjLayer):
				layer_indx = thisObjLayer[0]
				bpy.context.scene.layers[layer_indx] = True #this will set workspace at THIS layer
				break
		bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN') 
	#---------------#
	def invoke(self, context, event):
		self.func_unhide()
		self.func_unselect()
		pattern = str(bpy.context.scene.prop_obj_names)
		selected = self.func_select(pattern)
		bpy.types.Scene.data_obj_coords = [] #initliazed with each click		
		self.func_median_location(selected)

		for ob in selected:
			if ob.parent is None:
				bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(ob.location.x),str(ob.location.y),str(ob.location.z)])
			else:
				bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(ob.location.x),str(ob.location.y),str(ob.location.z)])

		
		return{"FINISHED"}
#--------------------------------------------------------------------------------
#									CLUSTERS
#--------------------------------------------------------------------------------
class OBJECTS_OT_calculate_clustering_param(bpy.types.Operator):
	bl_idname = "objects.calculate_clustering_param"
	bl_label = ""
	bl_description = "optomise DBSCAN results by calculating optimum values for step(2)"
	
	def invoke(self,context,event):
		gly_names = []
		gly_points = []
		
		bpy.types.Scene.data_sil_list = []
		bpy.types.Scene.data_nclusters = [] #initialise container with each new request
		bpy.types.Scene.data_epsx_list = []

		if not bpy.context.scene.data_obj_coords:
			populate_glycogen_coords()
		elif bpy.context.scene.prop_obj_names != "Glycogen" or bpy.context.scene.prop_obj_names != "Glycogen.":
			populate_glycogen_coords()
		
		for namePoint in bpy.context.scene.data_obj_coords:
			gly_names.append(namePoint[0])
			gly_points.append([float(point) for point in namePoint[2:]])

		np_points=np.array(gly_points)

		bpy.context.scene.prop_eps1= float("{0:.2f}".format(bpy.context.scene.prop_eps1))
		bpy.context.scene.prop_eps2= float("{0:.2f}".format(bpy.context.scene.prop_eps2))
		bpy.context.scene.prop_interval= float("{0:.2f}".format(bpy.context.scene.prop_interval))
		
		for x in np.arange(bpy.context.scene.prop_eps1, bpy.context.scene.prop_eps2, bpy.context.scene.prop_interval):
			db = DBSCAN(eps=x, min_samples=bpy.context.scene.prop_min_samples).fit(np_points)
			core_samples = db.core_sample_indices_

			labels = db.labels_
			no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
			sil_index = metrics.silhouette_score(np_points, labels)#sillohouette coeffecient			
			print('%0.2f 	%2d		%+0.3f' %(x, no_clusters,sil_index))
			
			bpy.types.Scene.data_sil_list.append(sil_index)
			bpy.types.Scene.data_nclusters.append(no_clusters)
			bpy.types.Scene.data_epsx_list.append(x)

		bpy.context.scene.prop_min_samples_s2 = bpy.context.scene.prop_min_samples
		bpy.context.scene.prop_optimum_eps = bpy.context.scene.data_epsx_list[bpy.context.scene.data_sil_list.index(max(bpy.context.scene.data_sil_list))]
		return{"FINISHED"}

class OBJECTS_OT_display_graph(bpy.types.Operator):
	bl_idname = "objects.display_graph"
	bl_label = ""
	bl_description = "dbscan result plot to graph"

	def invoke(self,context,event):
		#### NEEDS to start from a seperate process, to avoid blender from freezing
		import matplotlib
		from pylab import show
		from pylab import draw
		from pylab import plot
		from pylab import arange
		if bpy.types.Scene.data_sil_list:
			plot(arange(bpy.context.scene.prop_eps1, bpy.context.scene.prop_eps2, bpy.context.scene.prop_interval), bpy.types.Scene.data_sil_list)
			draw()
			show()
		#else error message
		return{"RUNNING_MODAL"}

class OBJECTS_OT_generate_clusters(bpy.types.Operator):
	bl_idname = "objects.generate_clusters"
	bl_label = ""
	bl_description = "step2: generate clusters + display in VIEW3D"

	objects_names = []
	np_points_ = []

	layer_index = 0
	layers_array = []

	def invoke(self,context,event):
		gly_names =[]
		gly_points = []
		bpy.types.Scene.data_glyc_clusters = []
		#in case user decided to skip optimum parameters calculations
		if not bpy.context.scene.data_obj_coords:
			populate_glycogen_coords()
		if bpy.context.scene.prop_obj_names != "Glycogen" or bpy.context.scene.prop_obj_names != "Glycogen.":
			populate_glycogen_coords()
		#DBSCAN
		for namePoint in bpy.context.scene.data_obj_coords:
			gly_names.append(namePoint[0])
			gly_points.append([float(point) for point in namePoint[2:]])
		np_points=np.array(gly_points)

		db = DBSCAN(eps=bpy.context.scene.prop_optimum_eps, min_samples=bpy.context.scene.prop_min_samples_s2).fit(np_points)
		core_samples = db.core_sample_indices_
		labels = db.labels_
		no_clusters = len(set(labels)) - (1 if -1 in labels else 0)
		sil_index = metrics.silhouette_score(np_points, labels)#sillohouette coeffecient			
		#print('%0.2f 	%2d		%+0.3f' %(x, no_clusters,sil_index))
		bpy.context.scene.prop_nclusters = str(len(set(labels)))
		bpy.context.scene.prop_silh = str(metrics.silhouette_score(np_points, labels))
		#Store clusters
		for lname, lpoints, labels in zip(gly_names,np_points,labels):
			bpy.types.Scene.data_glyc_clusters.append([(lname),([float(lpoints[0]),float(lpoints[1]),float(lpoints[2])]),(labels)])
		print('done generating clusters - next, drawing...')

							# --- Displaying Clusters in VIEW3D ---
		
		#initialise temp variables:
		grand_start = time.time()
		clusters_sorted_list = []
		cluster_points= []
		self.objects_names = []
		self.np_points_ = []
		_names = []
		#self.layer_index = 0
		#self.layers_array = []

		#1- sorting data_glyc_clusters
		clusters_sorted_list = sorted(bpy.context.scene.data_glyc_clusters ,key=lambda sortedCls:sortedCls[2])
		#[['Glycogen.176'], [3.1494250297546387, 2.828038215637207, 2.7060840129852295], 0]

		#loop for starting values and getting Glycogens layer number in the space:
		for row_data in clusters_sorted_list:
			if(row_data[2] != -1):
				self.get_glyco_layer(row_data[0])
				print(row_data[2])
				prev_label = row_data[2]
				break
		data_size = len(clusters_sorted_list)
		print(self.layer_index, self.layers_array)

		#generate a new colour material
		import random
		print(random.random())
		this_color = self.makeMaterial('color',
			(random.random(), random.random(), random.random())
			,(1,1,1),1)

		print("********************* START ***********************")
		start = time.time()
		for granule_info in clusters_sorted_list:
			data_size = data_size -1 
			if granule_info[2] != -1:
				dlabel = granule_info[2]
				dname = granule_info[0]
				_names.append(dname)
				dpoints = granule_info[1][0:]

				if prev_label == dlabel:
					#bpy.context.object.name = dname[0] WRONG!, instead:
					#variations in dname[0] or just dname
					objToSelect = bpy.data.objects[dname]
					objToSelect.select = True
					#on linux it was found that 'bpy.context.selected_objects[0]' might return another object which is NOT equal to objToSelect, therefore, we will loop:
					for correct_obj in bpy.context.selected_objects:
						if correct_obj == objToSelect:
							objToColor = correct_obj
							self.setMaterial(objToColor,this_color)
							break
					cluster_points.append([float(point) for point in dpoints[0:3]])
						
					self.objects_names.append(dname)
					self.np_points_ = np.array(cluster_points) #notgood overwrites in memory while loop continues
					if data_size == 0:
						end=time.time()
						print("read glycogen granules in:", end-start)
						ellip_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
						self.draw_ellipsoid(dlabel,ellip_color)
				else:
					end = time.time()
					print("read glycogen granules in:", end-start)
					ellip_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
					self.draw_ellipsoid(prev_label,ellip_color)
					
					start = time.time()
					this_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
					cluster_points = []
					self.objects_names = []
					prev_label = dlabel
					#bpy.context.object.name = dname[0] Wrong!, instead:
					objToSelect = bpy.data.objects[dname]
					objToSelect.select = True
					for correct_obj in bpy.context.selected_objects:
						if correct_obj == objToSelect:
							objToColor = correct_obj
							self.setMaterial(objToColor,this_color)
							break
					cluster_points.append([float(point) for point in dpoints[0:3]])
					self.objects_names.append(dname)
					self.np_points_ = np.array(cluster_points)

		#hiding other granules (NOISE), cleaning the worksapce view
		ellipsoidstr = "ellipsoid"
		#print(self.justnames)

		for noise_obj in bpy.data.objects:
			if noise_obj.name not in _names:
				noise_obj.select = False
				noise_obj.hide = True
				#print("hide", noise_obj.name)	
			else:
				#print("not hide", noise_obj.name)
				noise_obj.select = True
				noise_obj.hide = False

			match = re.search(ellipsoidstr+'*', noise_obj.name)
			if match:
				noise_obj.select = True
				noise_obj.hide = False
		
		print("Done! in: ", time.time()-grand_start)
		bpy.types.Scene.flag_clusters_measure = True
		bpy.types.Scene.flag_clusters_measured = False

		return{"FINISHED"}

	@classmethod
	## Get Glycogen Layer ##
	def get_glyco_layer(self, obj_name):
		obj = bpy.data.objects[obj_name]
		thisObjLayer = [i for i in range(len(obj.layers)) if obj.layers[i] == True] #returns objects layer
		if(thisObjLayer):
			self.layer_indx = thisObjLayer[0]
		
		for i in range(20): #0 to 19 total 20 layers
			if i == self.layer_indx:
				self.layers_array.append(True)
			else:
				self.layers_array.append(False)
		print(self.layer_indx)
		print(self.layers_array)
	###---- MAKE METERIAL FUNCTION----
	def makeMaterial(self, name, diffuse, specular, alpha):
		mat = bpy.data.materials.new(name)
		mat.diffuse_color = diffuse
		mat.diffuse_shader = 'LAMBERT'
		mat.diffuse_intensity = 1.0
		mat.specular_color = specular
		mat.specular_shader = 'COOKTORR'
		mat.specular_intensity = 0.5
		mat.alpha = alpha
		mat.ambient = 1
		return mat
	def setMaterial(self, ob, mat):
		#me = ob.data
		#me.materials.append(mat)
		ob.active_material = mat #-necessary for predrawn objects with pre active materials, otherwise, append is enough
		if ob.active_material.name == 'mat_1':
			print(ob.name, 'error: color not changed')
	###---- DRAW ELLIPSOID FUNCTION -----
	def draw_ellipsoid(self, dlabel, this_color):
		label = str(dlabel)
		tolerance = 0.01
		P = self.np_points_
		(N, d) = np.shape(P)
		d = float(d)        
		Q = np.vstack([np.copy(P.T), np.ones(N)]) 
		QT = Q.T
		# initializations
		err = 1.0 + tolerance
		u = (1.0 / N) * np.ones(N)
		# Khachiyan Algorithm
		while err > tolerance:
			V = np.dot(Q, np.dot(np.diag(u), QT))
			M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))   
			# M the diagonal vector of an NxN matrix
			j = np.argmax(M)
			maximum = M[j]
			step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
			new_u = (1.0 - step_size) * u
			new_u[j] += step_size
			err = np.linalg.norm(new_u - u)
			u = new_u
		#End of while                
		center = np.dot(P.T, u)
		# the A matrix for the ellipse
		A = linalg.inv(
					   np.dot(P.T, np.dot(np.diag(u), P)) - 
					   np.array([[a * b for b in center] for a in center])
					   ) / d                               
		U, s, rotation = linalg.svd(A)
		radii = 1.0/np.sqrt(s)
		u = np.linspace(0.0, 2.0 * np.pi, 100)
		v = np.linspace(0.0, np.pi, 100)                
		# cartesian coordinates that correspond to the spherical angles:
		x = radii[0] * np.outer(np.cos(u), np.sin(v))
		y = radii[1] * np.outer(np.sin(u), np.sin(v))
		z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
		# rotate accordingly
		start = time.time()
		for i in range(len(x)):
			for j in range(len(x)):
				[x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation)+center
		end = time.time()
		print("ellipsoid no.:",dlabel,"drawn in ", end-start, "no of granules:",len(self.objects_names))
		
		#we need a fake mesh (sphere)
		print(self.layers_array)
		bpy.ops.mesh.primitive_uv_sphere_add(size=0.03, location = (0.0,0.0,0.0),layers=self.layers_array)
		bpy.context.object.name = 'fake' + label
		#HERE last edit Sept10th,2014
		print(self.layer_indx)                                                                                                                            
		bpy.context.scene.layers[self.layer_indx] = True #11 but 12 in count
		obj = bpy.context.active_object
		if obj.name == 'fake'+label and obj.type == 'MESH':
			bpy.ops.object.mode_set(mode='EDIT')
			bpy.ops.mesh.delete(type='VERT') #delete all fake object's vertices
			me = obj.data
			bm = bmesh.from_edit_mesh(me)
			vertIndx=0
			#x,y,z start plotting vertices:                  
			for i in range(len(x)):
				for j in range(len(x)):
					bm.verts.new((x[i,j],y[i,j],z[i,j]))
			bpy.ops.object.mode_set(mode = 'OBJECT')        
			#using the shrinkWrap modifier,
			#size 3, to cover the size of any cluster. for the shrink-wrap to work properly as skin
			print(self.layers_array)
			bpy.ops.mesh.primitive_uv_sphere_add(size=3, location=center,layers=self.layers_array) 
			obj = bpy.context.active_object
			obj.name = 'ellipsoid'+label
			
			#adding ellipsoid name to objects_list so that it is set as parent to them
			self.objects_names.append(obj.name)
			
			self.setMaterial(bpy.context.object, this_color)        
			bpy.ops.object.modifier_add(type='SHRINKWRAP')
			bpy.context.object.modifiers['Shrinkwrap'].target = bpy.data.objects['fake'+label]
			bpy.context.object.modifiers['Shrinkwrap'].wrap_method = 'NEAREST_VERTEX'
			bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Shrinkwrap")
			bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='BOUNDS')

			#apply transperency to uv_sphere, then hide ico_sphere, give both different colours.
			obj.active_material.use_transparency = True
			obj.active_material.transparency_method='MASK'
			obj.active_material.alpha = 0.178571
			obj.show_transparent = True
			
			#make ellipsoid a parent
			self.makeParentOf()

			#to hide 'icosphere'
			bpy.ops.object.select_all(action='DESELECT')
			bpy.data.objects['fake'+label].hide = True
	###---- PARENT-CHILD FUNCTION -----
	def makeParentOf(self):
		children_list = self.objects_names
		counting = 0
		bpy.ops.object.select_all(action='DESELECT')
		obj_names = bpy.data.objects
		for child in children_list:
			for ob in obj_names:
				if ob.name == child:
					ob.select = True
					counting = counting + 1
					if counting == len(children_list):
						ob.select = True
						bpy.ops.object.parent_set(type='OBJECT')

	def init(self):
		bpy.types.Scene.data_glyc_clusters = []
		
class OBJECTS_OT_clusters_nearest_neighbours(bpy.types.Operator):
	bl_idname = "objects.clusters_nearest_neighbours"
	bl_label = ""
	bl_description = "output the center of the nearest neural object per cluster in 3D, This Button is activated by a flag"
	"""there is no filtering done in measurements. as most meshes wont include, pericytes or endothelial cells.
		the only filteration happens in excluding '-1' labeled clusters, and those are done automatically once
		clusters are generated.
	"""	
	def invoke(self,context,event):
		self.initialise()
		#1-compute clusters centroids:
		bpy.types.Scene.data_clusters_centroids = self.get_clusters_centroids()
		#2-compute distances between cluster centroid and rest of neural objects vertices
		#2-1 get bouton and spine centroids vertices:		
		patterns = []
		if "bouton" in bpy.context.scene.data_names:
			patterns.append('bouton')
		elif "Bouton" in bpy.context.scene.data_names:
			patterns.append('Bouton')
		if "spine" in bpy.context.scene.data_names:
			patterns.append('spine')
		elif "Spine" in bpy.context.scene.data_names:
			patterns.append('Spine')

		firstTime = True
		if patterns:
			for patt in patterns:
				temp_attrib_np, temp_verts_np = getVertices(patt, "Center Vertices")
				#print("getting Center vertices for: ", pattern)
				print("temp_attrib_np.shape for ",patt,temp_attrib_np.shape)
				print("temp_verts_np.shape for ",patt, temp_verts_np.shape)
				
				if firstTime:
					bpy.types.Scene.neur_obj_verts_np = temp_verts_np
					bpy.types.Scene.neur_obj_attrib_np = temp_attrib_np
					firstTime = False
					#to set the big array to same dimensions of concatinated arrays
				else:
					bpy.types.Scene.neur_obj_attrib_np =np.concatenate((bpy.types.Scene.neur_obj_attrib_np, temp_attrib_np),axis=0)
					bpy.types.Scene.neur_obj_verts_np = np.concatenate((bpy.types.Scene.neur_obj_verts_np, temp_verts_np),axis=0)
			
			#2-2 WE DONT get Endo and Pericyte All Vertices
			#2-3 extract only vertices from clusters centroids, then get the distances
			cls_centroids_vertices = []
			for i in bpy.types.Scene.data_clusters_centroids:
				cls_centroids_vertices.append([i[1],i[2],i[3]])
			cls_centroids_vertices = np.array(cls_centroids_vertices)
			
			result = get_closest_distance(cls_centroids_vertices,
				bpy.types.Scene.neur_obj_verts_np)
			
			for clusterId, objctsIndx, dist in result:
				bpy.types.Scene.data_clusters_distances.append(
					[np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,0],
					np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,1],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,0],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,1],
					dist
					])
				print(
					[np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,0],
					np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,1],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,0],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,1],
					dist
					])

		#if measurment array is populated correctly
		bpy.types.Scene.flag_clusters_measured = True

		return{"FINISHED"}
	@classmethod
	def initialise(self):
		bpy.types.Scene.data_clusters_centroids = []
		bpy.types.Scene.data_clusters_distances = []

	def get_clusters_centroids(self):
		glyNames = []
		glyPts = []
		clsLbls = []
		cls_centroids = []

		for glyName , points, clsLabel in bpy.context.scene.data_glyc_clusters:
			if clsLabel == -1:
				continue
			glyNames.append(glyName)
			clsLbls.append(clsLabel)
			glyPts.append([points[0],points[1],points[2]])
		#convert them to arrays
		clsLbls = np.array(clsLbls).astype(int)
		glyNames = np.array(glyNames)
		glyPts = np.array(glyPts).astype(float)

		occurences = itemfreq(clsLbls.astype(int))
		#occrences will give the total number of how many times each value occured 
		#e.g.,
		""" array([[ 0, 17],
       			  [ 1, 16]])  
		"""
		for i in occurences:
			indexes = np.where(clsLbls == i[0])
			centroid = np.mean(glyPts[indexes], axis=0)
			cls_centroids.append([i[0],centroid[0],centroid[1],centroid[2],i[1]])
			# cluster label, x,y,z, total#granules
		return cls_centroids
	
	"""def get_closest_distance(self,first_verts, second_verts):
		closests_points_info = []
		#first_verts = np.asarray(bpy.types.Scene.glycogen_verts_np)
		#second_verts = np.asarray(bpy.types.Scene.neur_obj_verts_np)
		first_verts = np.asarray(first_verts)
		second_verts = np.asarray(second_verts)
		
		#scipy cKDTree algorithm
		mytree = cKDTree(second_verts)
		dist, indexes = mytree.query(first_verts)
		dist= np.asarray(dist)
		indexes = np.asarray(indexes)
		a = np.vstack(dist)
		b = np.vstack(indexes)
		c = np.hstack((b,a))

		for i in range(0,len(c)):
			closests_points_info.append([i, c[i,0] , c[i,1]])
		
		return closests_points_info"""

class OBJECTS_OT_export_clusters_measures(bpy.types.Operator):
	bl_idname="objects.export_clusters_measures"
	bl_label= "write data to file"
	
	directory = bpy.props.StringProperty()
	filename = bpy.props.StringProperty()

	def execute(self,context):
		directory = self.directory
		filename = self.filename
		fname = filename+'.tab'
		the_name = os.path.join(directory, fname)
		
		if bpy.context.scene.prop_bool_clusters == True:
			print("writing data of [clusters centroids to spine-bouton centroids] distances")
			with open(the_name,'wt') as output:
				writer = csv.writer(output, delimiter='\t')
				writer.writerow(['CId','CCentroid','ObjectName','ObjectParent', 'Distance'])
				rowToWrite = []
				for CId, CCentroid,ObjectName,ObjectParent, Distance in bpy.types.Scene.data_clusters_distances:
					writer.writerow([CId, CCentroid,ObjectName,ObjectParent, Distance])
		
		elif bpy.context.scene.prop_bool_glycogen == True:
			print("writing data of [glycogen to closests neighbors objects- distances]")
			with open(the_name,'wt') as output:
				writer = csv.writer(output, delimiter='\t')
				writer.writerow(['Neural Object','No.Associated Glycogens','Glycogen Names','Distance'])

				rowToWrite = []
				for neurObj, glyList in bpy.types.Scene.neuro_glyList_dict.items():
					writeOnce = True
					for glyName in glyList:
						gdistance = [dist for gname, dist in bpy.types.Scene.data_glyc_to_neighb_dist if gname == glyName]
						if gdistance:
							if writeOnce:
								writer.writerow([neurObj, len(glyList), glyName, gdistance[0]])# gidstance[0]=0.091054856874, gdistance=[0.091054856874324699]
								writeOnce = False
							else:
								writer.writerow([ neurObj, 0 , glyName, gdistance[0]])
		return{"FINISHED"}
	def invoke(self,context,event):
		#print("bpy.context.scene.prop_bool_clusters",bpy.context.scene.prop_bool_clusters)
		WindowManager = context.window_manager
		WindowManager.fileselect_add(self)
		return{"RUNNING_MODAL"}

#--------------------------------------------------------------------------------
#                                   PANEL LAYOUT
#--------------------------------------------------------------------------------
class UI_VIEW3D_PT(bpy.types.Panel):
	bl_idname = "UI_VIEW3D_PT"
	bl_label = "Glycogen Analysis"
	bl_space_type = "VIEW_3D" #window type where the panel will be drawn, e.g.,"PROPERTIES"
	bl_region_type =  "UI"#"TOOL_PROPS" # "UI" , 
	#bl_context = "Object" 
	#bl_options = {"DEFAULT_CLOSED"}#http://www.blender.org/api/blender_python_api_2_70_5/bpy.types.Panel.html
	#bl_context left out, by default it wil be put on a new tab called "Misc" if bl_regions was UI, which scene_type panel should be showing
	#http://en.wikibooks.org/wiki/Blender_3D:_Noob_to_Pro/Advanced_Tutorials/Python_Scripting/Addon_User_Interface

	def draw(self,context):
		layout = self.layout
		scene = context.scene
		obj = context.object
		row1 = layout.row()
		row1.alignment = 'LEFT'
		row1.operator("view3d.activate_addon_button", "Activate/Reload addon")

		#hide-unhide objects -general panel box-
		box1 = layout.box()
		box1.label("---Hide/Unhide tool---")
		if bpy.context.scene.prop_obj_names:
			box1.prop(scene,'prop_obj_names')
			row1_box1 = box1.row(align=True)
			row1_box1.operator("view3d.display_selected", "Display")

		col1 = layout.column()
		col1.label("--- Measurements ---")
		row1_col1 = col1.row(align=True)
		#toggle two boolean variables (if they exist!)
		if bpy.types.Scene.prop_bool_glycogen:
			row1_col1.prop(scene,'prop_bool_glycogen')
			row1_col1.prop(scene,'prop_bool_clusters')

			if bpy.context.scene.prop_bool_clusters == True: 
				row2_layout = layout.row()
				row2_layout.alignment = 'LEFT'
				row2_layout.operator("VIEW3D_OT_show_info", icon = "INFO")

				if bpy.context.scene.prop_dbscan_info:
					infobox = layout.box()
					row2_layout= infobox.row()
					row2_layout.label("Clustering pareameters can either be set automaticlly")
					row2_info = infobox.row()
					row2_info.alignment = 'LEFT'
					row2_info.label("from step (1) or manually from step (2)")
				
				row3_step1 = layout.row(align=True)
				row3_step1.label("(1) Calculate optimum values for step(2):")
				row4_step1 = layout.row(align=True)
				row4_step1.prop(scene,'prop_min_samples')
				row5_step1 = layout.row(align=True)
				row5_step1.label("Epsilon Range:")
				row5_step1.prop(scene,'prop_eps1')
				row5_step1.prop(scene,'prop_eps2')
				row6_step1 = layout.row(align=True)
				row6_step1.prop(scene,'prop_interval')
				
				row7_step1 = layout.row()
				row7_step1.operator("objects.calculate_clustering_param", "Calculate")
				row7_step1.operator("objects.display_graph", "Display Graph")
				row7_step1.alignment = 'LEFT'
				
				#step2
				rowline = layout.row(align=True)
				rowline.label("----------------------------------------")
				row1_step2 = layout.row(align=True)
				row1_step2.label("(2) Generate Clusters:")
				
				row2_step2 = layout.row()
				row2_step2.prop(scene, 'prop_min_samples_s2')
				row3_step2 = layout.row()
				row3_step2.prop(scene, 'prop_optimum_eps')
				
				row4_step2 = layout.row()
				row4_step2.alignment = 'LEFT'
				row4_step2.operator("objects.generate_clusters", "Generate")

				row5_step2 = layout.row(align=True)
				row5_step2.label("---- Results --- :")
				
				row6_step2 = layout.row()
				row6_step2.enabled = False
				row6_step2.label("Number of Clusters:")
				row6_step2.prop(scene,"prop_nclusters")

				row7_step2 = layout.row()
				row7_step2.enabled = False
				row7_step2.label("Sillohouette Coefficient:")
				row7_step2.prop(scene,'prop_silh')
				
				row8_step2 = layout.row()
				row8_step2.operator("objects.clusters_nearest_neighbours", "Calculate Nearest Neighbour")
				row8_step2.enabled = False

				row9_step2 = layout.row()
				row9_step2.operator("objects.export_clusters_measures", "Export Measurements to File", icon = "DISK_DRIVE")
				row9_step2.enabled = False
				
				if bpy.context.scene.flag_clusters_measure:
					row8_step2.enabled = True

					if bpy.context.scene.flag_clusters_measured:
						row9_step2.enabled = True

			elif bpy.context.scene.prop_bool_glycogen == True:
				layout.operator("objects.glycogens_neareste_neighbors", "Calculate Nearest Neighbour")

				row_glyc = layout.row()
				row_glyc.operator("objects.export_clusters_measures", "Export Measurements to File", icon = "DISK_DRIVE")#export glycogen measures
				row_glyc.enabled = False
				
				if bpy.types.Scene.data_glyc_distances:
					row_glyc.enabled = True

					row1_glyc = layout.row()
					row1_glyc.label("Output:")

					row2_glyc = layout.row()
					row2_glyc.prop(scene,'prop_glyc_neighbours')

					row3_glyc = layout.row()
					row3_glyc.prop(scene,'prop_associated_glyco')

					row4_glyc = layout.row()
					row4_glyc.prop(scene,'prop_total_granules')
					row4_glyc.enabled = False
					row5_glyc = layout.row()
					row5_glyc.prop(scene,'prop_glyc_to_neighb_dist')
					row5_glyc.enabled = False
#--------------------------------------------------------------------------------
#                                   MAIN EXECUTE
#---------------------------------------------------------------------------------
def register():#register takes class name | panel takes bl_idname as string
	
	bpy.utils.register_module(__name__)
	#public variables
	bpy.types.Scene.data_names = None 
	bpy.types.Scene.prop_obj_names = None #dropdown list
	bpy.types.Scene.data_obj_coords = None
	bpy.types.Scene.prop_bool_glycogen = None #toggle checkbox
	bpy.types.Scene.prop_bool_clusters = None #toggle checkbox
	#dbscan info class
	bpy.types.Scene.prop_dbscan_info = False #flag
	#clustering option:
	bpy.types.Scene.prop_min_samples = None
	bpy.types.Scene.prop_eps1 = None
	bpy.types.Scene.prop_eps2 = None
	bpy.types.Scene.prop_interval = None
	bpy.types.Scene.data_sil_list = None
	#step2- generating clusters
	bpy.types.Scene.prop_min_samples_s2 = None #in sync with prop_min_samples step 1
	bpy.types.Scene.prop_optimum_eps = None
	bpy.types.Scene.data_glyc_clusters = None
	bpy.types.Scene.prop_nclusters = None
	bpy.types.Scene.prop_silh = None
	#clusters measurements:
	bpy.types.Scene.data_clusters_centroids = None
	bpy.types.Scene.flag_clusters_measure = None
	bpy.types.Scene.flag_clusters_measured = None
	bpy.types.Scene.data_clusters_distances = None
	#glycogens measurements:
	bpy.types.Scene.prop_glyc_neighbours = None
	bpy.types.Scene.prop_associated_glyco = None
	bpy.types.Scene.prop_total_granules = None
	bpy.types.Scene.prop_glyc_to_neighb_dist = None
	bpy.types.Scene.data_glyc_distances = None
	
def unregister():
	bpy.utils.unregister_module(__name__)
	#free memory
	if bpy.types.Scene.data_names: 
		del bpy.types.Scene.data_names
	if bpy.types.Scene.prop_obj_names: #dropdown list for storing objects names from the outline
		del bpy.types.Scene.prop_obj_names
	if bpy.types.Scene.data_obj_coords:
		del bpy.types.Scene.data_obj_coords
	if bpy.types.Scene.prop_bool_glycogen:
		del bpy.types.Scene.prop_bool_glycogen
	if bpy.types.Scene.prop_bool_clusters:
		del bpy.types.Scene.prop_bool_clusters
	#dbscan info class
	if bpy.types.Scene.prop_dbscan_info: #flag
		del bpy.types.Scene.prop_dbscan_info
	#clustering option:
	if bpy.types.Scene.prop_min_samples:
		del bpy.types.Scene.prop_min_samples
	if bpy.types.Scene.prop_eps1:
		del bpy.types.Scene.prop_eps1
	if bpy.types.Scene.prop_eps2:
		del bpy.types.Scene.prop_eps2
	if bpy.types.Scene.prop_interval:
		del bpy.types.Scene.prop_interval
	if bpy.types.Scene.data_sil_list:
		del bpy.types.Scene.data_sil_list
	if bpy.types.Scene.prop_min_samples_s2:
		del bpy.types.Scene.prop_min_samples_s2
	if bpy.types.Scene.prop_optimum_eps:
		del bpy.types.Scene.prop_optimum_eps
	if bpy.types.Scene.data_glyc_clusters:
		del bpy.types.Scene.data_glyc_clusters
	if bpy.types.Scene.prop_nclusters:
		del bpy.types.Scene.prop_nclusters
	if bpy.types.Scene.prop_silh:
		del bpy.types.Scene.prop_silh
	#clusters measurements:
	if bpy.types.Scene.data_clusters_centroids:
		del bpy.types.Scene.data_clusters_centroids
	if bpy.types.Scene.flag_clusters_measure:
		del bpy.types.Scene.flag_clusters_measure
	if bpy.types.Scene.flag_clusters_measured:
		del bpy.types.Scene.flag_clusters_measured
	if bpy.types.Scene.data_clusters_distances:
		del bpy.types.Scene.data_clusters_distances
	#glycogen measurements:
	if bpy.types.Scene.prop_glyc_neighbours:
		del bpy.types.Scene.prop_glyc_neighbours
	if bpy.types.Scene.prop_associated_glyco:
		del bpy.types.Scene.prop_associated_glyco
	if bpy.types.Scene.prop_total_granules:
		del bpy.types.Scene.prop_total_granules
	if bpy.types.Scene.prop_glyc_to_neighb_dist:
		del bpy.types.Scene.prop_glyc_to_neighb_dist
	if bpy.types.Scene.data_glyc_distances:
		del bpy.types.Scene.data_glyc_distances
	
	
	