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
import math
import mathutils

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
	bl_label = " Activate/Reset addon"

	def func_init(self):
		bpy.types.Scene.data_names = []
		bpy.types.Scene.prop_obj_names = bpy.props.EnumProperty(items=[])
		bpy.types.Scene.data_obj_coords = None
		bpy.types.Scene.data_glycogen_naming_convention = None # or 'glycogen' changes according to user
		
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
		if "Glycogen." in lst or "Glycogen" in lst:
			bpy.types.Scene.data_glycogen_naming_convention = "Glycogen"
		elif 'glycogen.' in lst or 'glycogen' in lst:
			bpy.types.Scene.data_glycogen_naming_convention = 'glycogen'
		
		if bpy.types.Scene.data_glycogen_naming_convention:
			print("Created")
			bpy.types.Scene.prop_bool_glycogen = bpy.props.BoolProperty(name="Glycogens", description=" ",
			update=update_prop_bool_glyc)
			bpy.types.Scene.prop_bool_clusters = bpy.props.BoolProperty(name="Clusters", description=" ",
			update=update_prop_bool_clust)
			bpy.types.Scene.prop_min_samples = bpy.props.IntProperty(name='Minimum Samples')
			bpy.types.Scene.prop_eps1 = bpy.props.FloatProperty(name='from')
			bpy.types.Scene.prop_eps2 = bpy.props.FloatProperty(name='to')
			bpy.types.Scene.prop_interval = bpy.props.FloatProperty(name='Interval')
			bpy.context.scene.prop_min_samples = 19
			bpy.context.scene.prop_eps1 = 0.35
			bpy.context.scene.prop_eps2 = 0.38
			bpy.context.scene.prop_interval = 0.01			
			#step 2
			bpy.types.Scene.prop_min_samples_s2 = bpy.props.IntProperty(name='Minimum Samples')
			bpy.types.Scene.prop_optimum_eps = bpy.props.FloatProperty(name='Optimum Epsilon')
			bpy.types.Scene.prop_nclusters = StringProperty(name='')
			bpy.types.Scene.prop_silh = StringProperty(name='')#Silhouette Coefficient
			#initialization:
			bpy.context.scene.prop_min_samples_s2 = 19
			bpy.context.scene.prop_optimum_eps = 0.0
			bpy.context.scene.prop_nclusters=""
			bpy.context.scene.prop_silh = ""
			#should initialise all data holders:
			bpy.types.Scene.data_glyc_clusters = [] #clustering-clusters section
			bpy.types.Scene.data_clusters_centroids = [] #clusters measures
			bpy.types.Scene.data_clusters_distances = [] # =
			bpy.types.Scene.flag_clusters_measure = False  # =
			bpy.types.Scene.flag_clusters_measured = False # = to activate export data button
			#glycogen Measures:
			empty = [("","","")]
			bpy.types.Scene.prop_glyc_neighbours = EnumProperty(name="Neibouring Objects",items=empty)
			bpy.types.Scene.prop_associated_glyco = EnumProperty(name='Associated Granules:',items=empty)
			bpy.types.Scene.prop_total_granules = StringProperty(name='Total Granules:')
			bpy.types.Scene.prop_glyc_to_neighb_dist = StringProperty(name='Distance:')
			bpy.types.Scene.data_glyc_distances = []

			bpy.types.Scene.data_glyc_distances_occur = []

			bpy.types.Scene.glycogen_attrib_np = np.array([])
			bpy.types.Scene.glycogen_verts_np = np.array([])
			bpy.types.Scene.neur_obj_verts_np = np.array([])
			bpy.types.Scene.neur_obj_attrib_np = np.array([])
			
			bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = {}
			bpy.types.Scene.neuro_glyList_dict = {}

			bpy.types.Scene.data_glyc_neighbours = []
			bpy.types.Scene.data_noOfGlyc = []
			bpy.types.Scene.data_glyc_to_neighb_dist = []

			#free memory from UIList gadeget:
			if bpy.types.Scene.UIList_glyc_neighb:
				index = len(bpy.context.scene.UIList_glyc_neighb)
				for i in range(len(bpy.context.scene.UIList_glyc_neighb)+1):
					list = bpy.context.scene.UIList_glyc_neighb
					list.remove(index)			

					if index >= 0:
						index = index - 1
					else:
						break
			
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
			bpy.types.Scene.prop_obj_names = bpy.props.EnumProperty(name='Objects Names',items=objects)
	#---------------#
	def invoke(self,context,event):
		self.func_init()
		self.func_updateObjectsNames()
		return{"FINISHED"}
#---------------------------------------------------------------------------------
# 	  					# UPDATE / PUBLIC FUNCTIONS #
#---------------------------------------------------------------------------------
def UIlistUpdate(item):
	listtemp = []
	#listtemp = bpy.types.Scene.neuro_glyList_dict[bpy.types.Scene.data_glyc_neighbours[scn.UIList_glyc_neighb_indx]]
	listtemp = bpy.types.Scene.neuro_glyList_dict[item.li_glyc_neighbours]
	enum2 = []
	for _index, enumval in enumerate(listtemp):
		enum2.append((enumval,enumval, ""))

	if len(enum2) == 5:
		strtemp = enum2[0]
		enum2.append(strtemp)

	#if item.li_glyc_neighbours != 'bouton1 Axon138_E':
	bpy.types.Scene.li_associated_glyco = EnumProperty(name="Associated Granules:", items=enum2)
	#else:
		#bpy.types.Scene.li_associated_glyco = EnumProperty(name="Associated Granules:", 
		#	items=[('Glycogen.1128', 'Glycogen.1128', ''), ('Glycogen.457', 'Glycogen.457', ''), 
		#	('Glycogen.448', 'Glycogen.448', ''), ('Glycogen.1119', 'Glycogen.1119', ''), 
		#	('Glycogen.1110', 'Glycogen.1110', '')])
	#	bpy.types.Scene.li_associated_glyco = EnumProperty(name="Associated Granules:", items=[('1','1','a'),
	#		('2','2','b'),('3','3','c'),('4','4','c'),('5','5','c')])
		
	#bpy.context.scene.li_associated_glyco = enum2[0][0]
	
	
	#update distance accordingly
	for glyname, dist in bpy.context.scene.data_glyc_to_neighb_dist:
		if glyname == bpy.context.scene.li_associated_glyco:
			bpy.context.scene.li_glyc_to_neighb_dist = str(dist)

	#return None

#def PTlistUpdate(self,context):
#	for glyname, dist in bpy.context.scene.data_glyc_to_neighb_dist:
#		if glyname == bpy.context.scene.li_associated_glyco:
#			bpy.context.scene.li_glyc_to_neighb_dist = str(dist)
#	return None

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
	# update3/3/15: adding granule dimension size
	func_unhide()
	func_unselect()
	selected = func_select(bpy.types.Scene.data_glycogen_naming_convention)
	print("populate_glycogen_coords: ",selected)
	bpy.types.Scene.data_obj_coords = []
	func_median_location(selected)

	for ob in selected:
		if ob.parent is None:
			bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(ob.location.x),str(ob.location.y),str(ob.location.z),str(ob.dimensions.x)])
		else:
			bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(ob.location.x),str(ob.location.y),str(ob.location.z),str(ob.dimensions.x)])

"""------ NEURO_MORPH --------"""
def cross_product(v0, v1):
    x =   v0[1]*v1[2] - v0[2]*v1[1]
    y = -(v0[0]*v1[2] - v0[2]*v1[0])
    z =   v0[0]*v1[1] - v0[1]*v1[0]
    return [x,y,z]
def dot_product(v0,v1):
    vec = [v0[n]*v1[n] for n in range(len(v0))]
    return sum(vec)
# get area of single triangle
def get_area_tri(tri):
    # tri = [p0, p1, p2]
    p0 = tri[0]
    p1 = tri[1]
    p2 = tri[2]
    area = mathutils.geometry.area_tri(p0, p1, p2)
    return area
# get signed volume contribution from single triangle
def get_vol_tri(tri):  
    # tri = [p0, p1, p2],  pn = [x, y, z]
    p0 = tri[0]
    p1 = tri[1]
    p2 = tri[2]
    vcross = cross_product(p1,p2)
    vdot = dot_product(p0, vcross)
    vol = vdot/6
    return vol

def fget_SA(this_ob):
    obj = this_ob.data
    if hasattr(obj, 'polygons'):
        n_faces = len(obj.polygons)
        SA = 0
        for f in range(0, n_faces):
            n_vertices = len(obj.polygons[f].vertices[:])
            if n_vertices != 3:  # faces must be triangles
                return 'highlight a subregion'
            tri = [0] * n_vertices
            for v in range(0, n_vertices):
                tri[v] = obj.vertices[obj.polygons[f].vertices[v]].co
            SA += get_area_tri(tri)
        return SA
    else:
        return 'property not available'
def fget_vol(this_ob):
    obj = this_ob.data
    if hasattr(obj, 'polygons'):
    	# dani: this part is commented out, as we allready have solid objects that are meant to be measured for volume
        # if mesh not closed, don't calculate volume
        #if this_ob.is_open:
        #    return 'open mesh has no volume'
        #else:
        n_faces = len(obj.polygons)
        vol = 0
        for f in range(0, n_faces):
            n_vertices = len(obj.polygons[f].vertices[:])
            if n_vertices != 3:  # faces must be triangles
                for v in range(0, n_vertices):
                    return 'highlight a subregion'
            tri = [0] * n_vertices
            for v in range(0, n_vertices):
                tri[v] = obj.vertices[obj.polygons[f].vertices[v]].co
            vol += get_vol_tri(tri)
        return vol
    else:
        return 'property not available'
# these functions need to exist but are not used
def fset_vol(self, value):
    self.vol = -1
def fset_SA(self, value):
    self.SA = -1
"""----------end--------------"""
def paths():
    # RELEASE SCRIPTS: official scripts distributed in Blender releases
    addon_paths = bpy.utils.script_paths("addons")

    # CONTRIB SCRIPTS: good for testing but not official scripts yet
    # if folder addons_contrib/ exists, scripts in there will be loaded too
    addon_paths += bpy.utils.script_paths("addons_contrib")

    # EXTERN SCRIPTS: external projects scripts
    # if folder addons_extern/ exists, scripts in there will be loaded too
    addon_paths += bpy.utils.script_paths("addons_extern")

    return addon_paths


def glyc_toGlobal_coords(ob_index,objs_attrib,objs_verts,objs_surf,objs_solids,kwords_flg):
	# for boutons and spines case: surf_areas will be calculated from surfs, volume from solids or vols
	# we check if neuro_morph is uploaded (recommended), otherwise, we compute the two values 
	path_list = paths()

	if kwords_flg: # special case
		for path in path_list:
			for mod_name, mod_path in bpy.path.module_names(path):
				if mod_name == 'NeuroMorph Measurement Tools:  Submesh Volume, Surface Area, and Length':
					surf_are_time_start = time.time()
					# SA and VOl computed from neuroMorph
					bpy.tyeps.Scene.surf_area = objs_surf[ob_index].SA
					if objs_solids[ob_index] != 'incomplete_cell':
						bpy.types.Scene.volume = objs_solids[ob_index].vol
					else:
						bpy.types.Scene.volume = 'incomplete_cell' #0
					break
				else:
					surf_are_time_start = time.time()
					bpy.types.Scene.surf_area = fget_SA(objs_surf[ob_index])#property(fget_SA, fset_SA)
					if objs_solids[ob_index] != 'incomplete_cell':
						bpy.types.Scene.volume = fget_vol(objs_solids[ob_index])
					else:
						bpy.types.Scene.volume = 'incomplete_cell' #0
					break
		
		coord_matrix = objs_surf[ob_index].matrix_world
		
		# this loop and saves all vertices for one object
		for v in objs_surf[ob_index].data.vertices:
			loc = v.co
			globalCoords = coord_matrix *loc
			#object has no parent
			if objs_surf[ob_index].parent is None:
				objs_verts.append([str(globalCoords.x)
					,str(globalCoords.y)
					,str(globalCoords.z) ])
				# object is a bouton/spine
				if kwords_flg:
					objs_attrib.append([objs_surf[ob_index].name
						,'None'
						, bpy.types.Scene.surf_area
						, bpy.types.Scene.volume ])
				# object is endo or pericyte
				else:
					obj = ob_index #for non bouton/spine objects, the ob_index is actually an object. so we correct variable name to obj
					objs_attrib.append([obj.name,'None'])
			# object has a parent
			else:
				objs_verts.append([str(globalCoords.x)
						,str(globalCoords.y)
						,str(globalCoords.z) ])
				if kwords_flg:
					objs_attrib.append([objs_surf[ob_index].name
						,objs_surf[ob_index].parent.name
						,bpy.types.Scene.surf_area
						,bpy.types.Scene.volume ])
				else:
					obj = ob_index
					objs_attrib.append([obj.name, obj.parent.name])
	else: # Default case. e.g., for Endo's and pericytes
		obj = ob_index
		bpy.types.Scene.surf_area = 'None' # we need this to match the size of the numpy array
		bpy.types.Scene.volume = 'None'

		coord_matrix = obj.matrix_world
		# this loop and saves all vertices for one object
		for v in obj.data.vertices:
			loc = v.co
			globalCoords = coord_matrix *loc
			#object has no parent
			if obj.parent is None:
				objs_verts.append([str(globalCoords.x)
					,str(globalCoords.y)
					,str(globalCoords.z) ])
				objs_attrib.append([obj.name
					,'None'
					, bpy.types.Scene.surf_area
					, bpy.types.Scene.volume ])
			# object has a parent
			else:
				objs_verts.append([str(globalCoords.x)
						,str(globalCoords.y)
						,str(globalCoords.z) ])
				objs_attrib.append([obj.name
					,obj.parent.name
					,bpy.types.Scene.surf_area
					,bpy.types.Scene.volume ])
	return(objs_attrib,objs_verts)

# getVertices will seperate vertices from attributes
def getVertices(pattern,coords_type): 
	func_unhide()
	func_unselect()
	# These lists used to store a whole Objects
	selected_objects = [] #glycogens and non neuronal objects (e.g., Endo's and Pericytes)
	selected_objects_solids=[] # solids and surfs lists are for boutons and spines
	selected_objects_surfs=[] # Those 3 should all be with the same length, stores boutons and spines
	
	objs_attrib = [] #either name, parent or just name
	objs_verts = []
	objs_attrib_matrix = []
	objs_verts_matrix = []
	kwords = ['bouton','spine','Bouton','Spine']
	kwords_flg = False

	if pattern in kwords:
		kwords_flg = True #This category has an area and volume (for spines and boutons)
		is_incomplete = False
	
	#================ storing objects according to their type ===================
	for ob in bpy.data.objects:
		if ob.type != 'MESH': #only taking MESHS, no CURVES,..etc
			continue
		
		if kwords_flg: #flag for a bouton or spine object
			# an incomplete object is detected either by finding a *vol | *solid copy of it. or from its object.name(), where it starts by "in"
			if re.search(pattern+'.*surf*', ob.name): #for every surf there's either a solid or vol object depends on how they're named from neuroMorph
				is_incomplete = True # an incomplete object will only have surf, unless found otherwise, by finding a solid or volume
				ob.select=True
				selected_objects_surfs.append(ob)
				# search in children list of that object's parent. to extract solids or volume objects
				one_word=ob.name.rsplit('_',1)

				for child in ob.parent.children:
					# there can either be a solid or vol per parent. not correct to have both: wrong! in case if there was 3 children, one has a solid (complete) and one doesnt (incomplete)
					# update Mar9th-15: solids exist if an object was complete. 
					# for incomplete objects there's only a surface area value no volumes

					this_child=child.name.rsplit('_',1) # in case if a parent has one complete and another incomplete, hence 3 children, 
														# we need to match first part of object_name to recognise its copies
					if one_word[0] == this_child[0]:
						if re.search(one_word[0]+ '.*solid*' + '|' + one_word[0] + '.*volu*', child.name):
							is_incomplete = False  #object is complete
							selected_objects_solids.append(child)
							break
				if is_incomplete:
					selected_objects_solids.append('incomplete_cell')					
			else:
				ob.hide=True
				continue
		else:
			# neither bouton nor spine
			if re.search(pattern+'*', ob.name):
				ob.select = True
				selected_objects.append(ob)
			else:
				ob.hide=True
				continue
	#====================================================
	
	if selected_objects:
		print("selected_objects",len(selected_objects), selected_objects[0])
	
	# Total SOLID/VOLUME objects = total SURF objects
	if selected_objects_solids:
		print("selected_objects_solids",len(selected_objects_solids),selected_objects_solids[0])
	if selected_objects_surfs and selected_objects_solids:
		print("selected_objects_surfs",len(selected_objects_surfs),selected_objects_surfs[0])
	
	#================ vertices/attributes ===================

	if bpy.context.scene.prop_bool_glycogen:
		
		if selected_objects and coords_type == "Centroids": # usual case for glycogens
			func_median_location(selected_objects)
			
			for ob in selected_objects:
				objs_verts.append([str(ob.location.x),str(ob.location.y),str(ob.location.z)])
				
				if ob.parent is None:
					objs_attrib.append([ob.name, 'None', ob.dimensions.x]) # we need to get size dimensions per granule
				else:
					objs_attrib.append([ob.name, ob.parent.name , ob.dimensions.x])
				
		elif coords_type == 'All Vertices':
			# usual case for boutons and spines
			if kwords_flg and selected_objects_surfs:
				time_for_getting_solid_objects=time.time()
				# we take the index
				for ob_index in range(len(selected_objects_surfs)):
					objs_attrib, objs_verts = glyc_toGlobal_coords(ob_index
						,objs_attrib
						,objs_verts
						,selected_objects_surfs
						,selected_objects_solids
						,kwords_flg)
				print("getting solids/vols for ",pattern, 'Done in:', time.time()-time_for_getting_solid_objects)
			# case for Endo's and Pericytes
			elif not kwords_flg and selected_objects:
				for ob in selected_objects:
					objs_attrib, objs_verts = glyc_toGlobal_coords(ob
						,objs_attrib
						,objs_verts
						,[]
						,[]
						,kwords_flg)
		#========================================================

	elif  bpy.context.scene.prop_bool_clusters:
		if selected_objects_surfs and selected_objects_solids:# and coords_type == "Centroids":
			# usual case for spines and boutons
			#if kwords_flg and selected_objects_surfs:
			if coords_type == "Centroids":
			
				is_neuro_morpth = False
				path_list = paths()

				# ----------- checking for neuro morph -------
				for path in path_list:
					for mod_name, mod_path in bpy.path.module_names(path):
						if mod_name == 'NeuroMorph Measurement Tools:  Submesh Volume, Surface Area, and Length':
							is_neuro_morpth = True
							break
				# --------------------------------------------
				func_median_location(selected_objects_surfs)

				the_index = 0
				for ob in selected_objects_surfs:
					objs_verts.append([str(ob.location.x)
						,str(ob.location.y)
						,str(ob.location.z)])
	
					if is_neuro_morpth:
						bpy.types.Scene.surf_area = ob.SA

						if selected_objects_solids[the_index] == 'incomplete_cell':						
							bpy.types.Scene.volume = 'incomplete_cell'
						else:
							bpy.types.Scene.volume = ob.vol
					else:
						bpy.types.Scene.surf_area = fget_SA(ob)
						if selected_objects_solids[the_index] == 'incomplete_cell':
							bpy.types.Scene.volume = 'incomplete_cell'
						else:
							bpy.types.Scene.volume = fget_vol(ob)	
					
					if ob.parent is None:
						objs_attrib.append([ob.name
							,'None'
							, bpy.types.Scene.surf_area
							, bpy.types.Scene.volume ])
					else:
						objs_attrib.append([ob.name
							, ob.parent.name
							, bpy.types.Scene.surf_area
							, bpy.types.Scene.volume ])
					the_index = the_index + 1
			elif coords_type == 'All Vertices':
				#case of spines & boutons All vertices
				time_for_getting_solid_objects=time.time()
				# we take the index
				for ob_index in range(len(selected_objects_surfs)):
					objs_attrib, objs_verts = glyc_toGlobal_coords(ob_index
						,objs_attrib
						,objs_verts
						,selected_objects_surfs
						,selected_objects_solids
						,kwords_flg)
				print("getting solids/vols for ",pattern, 'Done in:', time.time()-time_for_getting_solid_objects)
		# case for Endo's and Pericytes
		elif not kwords_flg and selected_objects:
			if coords_type == 'Centroids':
				func_median_location(selected_objects)
				
				for ob in selected_objects:
					objs_verts.append([str(ob.location.x)
						,str(ob.location.y)
						,str(ob.location.z)])
					if ob.parent is None:
						objs_attrib.append([ob.name
							,'None'
							,'None_sa'
							,'None_vol'])
					else:
						objs_attrib.append([ob.name
							,ob.parent.name
							,'None_sa'
							,'None_vol'])

			elif coords_type == 'All Vertices':
				for ob in selected_objects:
					objs_attrib, objs_verts = glyc_toGlobal_coords(ob
						,objs_attrib
						,objs_verts
						,[]
						,[]
						,kwords_flg)

	return(np.array(objs_attrib),np.array(objs_verts))

def get_solid_objects(childWparentList):
	list_result = []
	printout=[]
	found = False
	for child, parent in childWparentList:
		found = False
		n = child.rsplit('_',1)#take the frist syllabes. e.g., bouton1_surf.059. we take bouton1_ only
		
		for this_child in bpy.data.objects[parent].children:
			match=re.search(n[0]+'.solid*'+'|'+ n[0]+'.volu*', this_child.name)# '*' 0 or more repetition
			if match:
				#print('this_child, parent',this_child.name, parent)
				list_result.append([this_child.name, parent])
				found = True
				break

		if not found:
			printout.append(parent)
	
	print4check(printout)
	return list_result

def print4check(lista):
	templist=[]
	countInstances = Counter(lista) #countInstances is a dictionary
	#{'listItem': itemfreq()...}
	for neur_obj_name, instances in countInstances.items():
		templist.append([neur_obj_name,instances])
	print('missed solids:',templist)

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
	#[index(glycogen#),objectindex,distance]?
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
		bpy.types.Scene.glycogen_attrib_np, bpy.types.Scene.glycogen_verts_np = getVertices(
			bpy.types.Scene.data_glycogen_naming_convention, "Centroids")

		patterns = [] # we can switch to a dictionary of options
		if "bouton" in bpy.context.scene.data_names:
			patterns.append('bouton')
		elif "Bouton" in bpy.context.scene.data_names:
			patterns.append('Bouton')
		if "spine" in bpy.context.scene.data_names:
			patterns.append('spine')
		elif "Spine" in bpy.context.scene.data_names:
			patterns.append('Spine')
		if 'Endothelial cell ' in bpy.context.scene.data_names:
			patterns.append('Endothelial cell')
		elif 'endothelial cell' in bpy.context.scene.data_names:
			patterns.append('endothelial cell')
		if 'Pericyte' in bpy.context.scene.data_names:
			patterns.append('Pericyte')
		elif 'pericyte' in bpy.context.scene.data_names:
			patterns.append('pericyte')

		firstTime = True

		for patt in patterns:
			print("getting vertices for: ", patt)
			temp_attrib_np, temp_verts_np = getVertices(patt, "All Vertices")
			if temp_attrib_np.size:
				if firstTime:
					bpy.types.Scene.neur_obj_verts_np = temp_verts_np
					bpy.types.Scene.neur_obj_attrib_np = temp_attrib_np
					firstTime = False
					#to set the big array to same dimensions of concatinated arrays
				else:
					bpy.types.Scene.neur_obj_attrib_np =np.concatenate((bpy.types.Scene.neur_obj_attrib_np, temp_attrib_np),axis=0)
					bpy.types.Scene.neur_obj_verts_np = np.concatenate((bpy.types.Scene.neur_obj_verts_np, temp_verts_np),axis=0)

		print("this attrib: size ",patt,bpy.types.Scene.neur_obj_attrib_np.size)
		print("this attrib: len ",patt,len(bpy.types.Scene.neur_obj_attrib_np))
		print("this vertices: size ",patt,bpy.types.Scene.neur_obj_verts_np.size)
		print("this vertices: len ",patt,len(bpy.types.Scene.neur_obj_verts_np))

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
		#-5 UIList populate
		countIndx = 0
		for neural_obj_name, noOfGlycogens in bpy.types.Scene.data_glyc_distances_occur:
			my_item = bpy.context.scene.UIList_glyc_neighb.add()
			my_item.li_glyc_neighbours=neural_obj_name
			my_item.li_total_granules = str(noOfGlycogens)
			countIndx+=1
			my_item.li_row_number = str(countIndx)

	def initialise(self):
		
		bpy.types.Scene.data_glyc_distances = []
		bpy.types.Scene.data_glyc_distances_occur = []
		bpy.types.Scene.li_associated_glyco = EnumProperty(name='Associated Granules:',items=[])
		bpy.types.Scene.li_glyc_to_neighb_dist = StringProperty(name='Distance:',default="")

	def occurences(self):
		objects_names=[]
		objects_SA_vol=[]
		gly_names=[]
		gly_sizes=[]

		bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = {} # shoud be used for data export
		bpy.types.Scene.neuro_glyList_dict = {} #used for UI display (a group by done on 'sorted')
		dict_temp1 = {} # needed to sort data by object name, result in bpy.types.Scene.neuro_gly_glyFreq_dic_sorted
		bpy.types.Scene.dict_temp2= {} # stores glynames with sizes {key:glyname,val:glysize}
		bpy.types.Scene.dict_temp3= {} # stores objectsNames with sa and vol {key:objName,val:<sa, vol>}
		
		### +++++ all datastructures indexes correspond with the main glycogen index++++++ ####
		
		closest_points_np = (np.array(bpy.types.Scene.data_glyc_distances)).astype(int)
		for k in range(0, len(closest_points_np)):
			#store closests neighbours (parent  child)
			objects_names.append(" ".join
				((
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],0],
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],1]
				))
				)#join (object name, parent) with a " " between them

			#store size for each glycogen
			gly_names.append(bpy.types.Scene.glycogen_attrib_np[closest_points_np[k,0],0])
			gly_sizes.append(bpy.types.Scene.glycogen_attrib_np[closest_points_np[k,0],2]) #indx 1 is the parent
			
		#this object is either a spine or bouton, has a surface area and volume attributes
		#store volume and surface area
			if bpy.types.Scene.neur_obj_attrib_np.shape[1] == 4:
				objects_SA_vol.append(" ".join
					((
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],2],
					bpy.types.Scene.neur_obj_attrib_np[closest_points_np[k,1],3]
					))
					)
		print('len(objects_names)',len(objects_names))
		print('len(gly_names)',len(gly_names))
		print('len(objects_SA_vol', len(objects_SA_vol))

		if not objects_names or not gly_names or not objects_SA_vol:
			print('error, faulty array - occurences function- closest neighbours glycoges')

		'''now we can create dictionary from obj&gly names:'''
		''' one nested dictionary () '''
		for i in range(0,len(gly_names)):
			#the below method in populating dictionaries will perform (itemfreq) on objects names so that it will not duplicate keys (robust?)
			dict_temp1[gly_names[i]] = objects_names[i] #we only need to sort this
			bpy.types.Scene.dict_temp2[gly_names[i]] = gly_sizes[i] #public
			
			#we need {objects_names, objects_SA_Vol} without repetition, its len < len(glycogens). wrong to place it in a glyocgen loop
			bpy.types.Scene.dict_temp3[objects_names[i]] = objects_SA_vol[i] #public
		
		#sort the dictionary:b = OrderedDict(sorted(a.items()))
		bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = OrderedDict(sorted(dict_temp1.items(), key=lambda val: val[1])) # sort by neural object name

		#now we need to group glycogens in a list per neuro
		#===============================================
		#Populating bpy.types.Scene.neuro_glyList_dict
		templist1 = []
		current_neuro = list(bpy.types.Scene.neuro_gly_glyFreq_dic_sorted.values())[0]
		for glyname, objname in bpy.types.Scene.neuro_gly_glyFreq_dic_sorted.items(): # .item refers to a pair(key,value), switching keys to values and values to keys
			if objname == str(current_neuro):
				templist1.append(glyname)
			else:
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

		return templist


#=========+++++++ UIList Gadget +++++++===========+
# Custom properties, will be saved with .blend file.
class PG_List_Entry(bpy.types.PropertyGroup):
	li_glyc_neighbours = StringProperty(name="Object Name")
	li_total_granules = StringProperty(name='Total Granules')
	li_row_number = StringProperty(name='')

class SCENE_UI_List(bpy.types.UIList):
	def draw_item(self,context,layout, data, item, icon, active_data, active_propname, index, flt_flag):
		custom_icon = 'OBJECT_DATAMODE'
		if self.layout_type in {'DEFAULT','COMPACT'}:
			#col = layout.column()
			#col.alignment = 'LEFT'
			#col.prop(item,'li_row_number')
			#col.enabled = False
			layout.label(item.li_glyc_neighbours)
		elif self.layout_type in {'GRID'}:
			layout.alignment = 'CENTER'
			layout.label("", icon=custom_icon)
	def filter_items(self, context, data, propname):
		#print('self.filter_name',self.filter_name)
		pgroup = getattr(data,propname)
		helper_funcs=bpy.types.UI_UL_list
		# default return values
		flt_flags = []
		flt_neworder=[]
		# Filtering by name
		if self.filter_name:
			flt_flags = helper_funcs.filter_items_by_name(self.filter_name, self.bitflag_filter_item, pgroup, "li_glyc_neighbours")
		if not flt_flags:
			flt_flags = [self.bitflag_filter_item] * len(pgroup)
			flt_neworder = helper_funcs.sort_items_by_name(pgroup, "li_glyc_neighbours")
		#print('flt_flags',flt_flags)
		#print('flt_neworder',flt_neworder)
		return flt_flags, flt_neworder
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
			if len(pattern) == 1: # for naming conventions like synapses= d0s0a143b1_E
				match = re.search('\A'+pattern+'\d',ob.name)
			else:
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
		# update3/3/15: adding granule dimension size
		self.func_unhide()
		self.func_unselect()
		pattern = str(bpy.context.scene.prop_obj_names)
		selected = self.func_select(pattern)
		bpy.types.Scene.data_obj_coords = [] #initliazed with each click		
		self.func_median_location(selected)

		for ob in selected:
			if ob.parent is None:
				if bpy.context.scene.prop_obj_names == bpy.types.Scene.data_glycogen_naming_convention:# "Glycogen":
					bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(ob.location.x),str(ob.location.y),str(ob.location.z),str(ob.dimensions.x)])
				else:
					bpy.types.Scene.data_obj_coords.append([ob.name, 'none', str(ob.location.x),str(ob.location.y),str(ob.location.z)])
			else:
				if bpy.context.scene.prop_obj_names == bpy.types.Scene.data_glycogen_naming_convention:# "Glycogen":
					bpy.types.Scene.data_obj_coords.append([ob.name, ob.parent.name , str(ob.location.x),str(ob.location.y),str(ob.location.z), str(ob.dimensions.x)])
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

		if not bpy.context.scene.data_obj_coords:#if its empty, fill it will glycogens info
			populate_glycogen_coords()
		elif bpy.context.scene.prop_obj_names != bpy.types.Scene.data_glycogen_naming_convention or bpy.context.scene.prop_obj_names != bpy.types.Scene.data_glycogen_naming_convention +'.':
		#elif bpy.context.scene.prop_obj_names != "Glycogen" or bpy.context.scene.prop_obj_names != "Glycogen.":
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
		# necessary init
		gly_names =[]
		gly_points = []
		bpy.types.Scene.data_glyc_clusters = []
		
		#in case user decided to skip optimum parameters calculations
		if not bpy.context.scene.data_obj_coords:
			populate_glycogen_coords()
		if bpy.context.scene.prop_obj_names != bpy.types.Scene.data_glycogen_naming_convention or bpy.context.scene.prop_obj_names != bpy.types.Scene.data_glycogen_naming_convention +'.':
		#if bpy.context.scene.prop_obj_names != "Glycogen" or bpy.context.scene.prop_obj_names != "Glycogen.":
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
				print(row_data[2])#should be zero
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
					objToSelect = bpy.data.objects[dname]
					objToSelect.select = True
					for correct_obj in bpy.context.selected_objects:
						if correct_obj == objToSelect:
							objToColor = correct_obj
							self.setMaterial(objToColor,this_color)
							break
					# jun2nd-15:
					pointsList=[]
					i=1
					for point in dpoints:
						dexp=len(str(point).split('.')[1])
						if i==3:#change zcoords only
							pointsList.append('{0:.{precision}f}'.format(point,precision=dexp-1))
						else:
							pointsList.append('{0:.{precision}f}'.format(point,precision=dexp))
						i=i+1
					cluster_points.append(pointsList)
					#cluster_points.append([float(point) for point in dpoints[0:3]])
					self.np_points_ = np.array(cluster_points, dtype=np.float64)
					# end jun2nd changes

					self.objects_names.append(dname)

					if data_size == 0:
						end=time.time()
						print("read glycogen granules in:", end-start)
						ellip_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
						#jun4th
						self.call_draw_ellipsoid(dlabel, ellip_color)
						#self.draw_ellipsoid(dlabel,ellip_color)
						#end jun4th
				else:
					end = time.time()
					print("read glycogen granules in:", end-start)
					ellip_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
					# jun2nd-15
					self.call_draw_ellipsoid(prev_label, ellip_color)
					#self.draw_ellipsoid(prev_label,ellip_color)
					# end
					start = time.time()
					this_color = self.makeMaterial('color', (random.random(),random.random(),random.random()),(1,1,1), 1)
					cluster_points = []
					self.objects_names = []
					prev_label = dlabel
					objToSelect = bpy.data.objects[dname]
					objToSelect.select = True
					for correct_obj in bpy.context.selected_objects:
						if correct_obj == objToSelect:
							objToColor = correct_obj
							self.setMaterial(objToColor,this_color)
							break
					# Jun2nd-15:
					pointsList=[]
					i=1
					for point in dpoints:
						dexp=len(str(point).split('.')[1])
						if i==3:# we want zCoords only
							pointsList.append('{0:.{precision}f}'.format(point,precision=dexp-1))
						else:
							pointsList.append('{0:.{precision}f}'.format(point,precision=dexp))
						i=i+1
					cluster_points.append(pointsList)
					self.np_points_ = np.array(cluster_points,dtype=np.float64)
					#cluster_points.append([float(point) for point in dpoints[0:3]])
					# end jun2nd change
					
					self.objects_names.append(dname)
					
		#hiding other granules (NOISE), cleaning the worksapce view
		ellipsoidstr = "ellipsoid"
		
		for noise_obj in bpy.data.objects:
			if noise_obj.name not in _names:
				noise_obj.select = False
				noise_obj.hide = True
			else:
				noise_obj.select = True
				noise_obj.hide = False

			match = re.search(ellipsoidstr+'*', noise_obj.name)
			if match:
				noise_obj.select = True
				noise_obj.hide = False
		
		print("Done! in: ", time.time()-grand_start)
		bpy.types.Scene.flag_clusters_measure = True
		bpy.types.Scene.flag_clusters_measured = False

		# hiding relationship lines if they're ON
		self.relationship_lines_off()
		# self.relationship_lines_off doesnt work without bracekts in python3, doesnt give an error either!

		return{"FINISHED"}

	## Get Glycogen Layer ##
	@classmethod
	def get_glyco_layer(self, obj_name):
		self.layer_index = 0
		self.layers_array = []
		obj = bpy.data.objects[obj_name]
		thisObjLayer = [i for i in range(len(obj.layers)) if obj.layers[i] == True] #returns objects layer
		if(thisObjLayer):
			self.layer_index = thisObjLayer[0]
		
		for i in range(20): #0 to 19 total 20 layers
			if i == self.layer_index:
				self.layers_array.append(True)
			else:
				self.layers_array.append(False)
		print(self.layer_index)
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
	def call_draw_ellipsoid(self,dlabel, this_color):
		# error_flag to avoid singular matrix inverse calculations April29th
		error_flag=True
		fix_flag=False
		err_rounds = 0
		# extract number of decimal places to add noise to the right most place of the first value in the points matrix.
		# this is to avoid singular matrixes from popping up
		dexp = len(str(self.np_points_[0,2]).split('.')[1])
		noise = 1/(10**dexp)

		while error_flag:
			if fix_flag:
				err_rounds = err_rounds + 1
				old_value = self.np_points_[0,2]
				self.np_points_[0,2] = self.np_points_[0,2] + noise #least amount of noise
				print('for ellipsoid ', dlabel, 'fix_flag wass applied for the ',err_rounds,'time. Matrix element was ',
					old_value ,' changed to ',self.np_points_[0,2],'noise value approximte: ',noise)
				fix_flag = False #switch it off after applied fix

			error_flag, fix_flag = self.draw_ellipsoid(dlabel, this_color)
			if not error_flag:
				break
			print('for ellipsoid ', dlabel, 'draw_ellipsoid has been called for the ',err_rounds,
				'time. Matrix element ', self.np_points_[0,2],'noise value is: ',noise)

	def draw_ellipsoid(self, dlabel, this_color):
		label = str(dlabel)

			# initializations
		tolerance = 0.01
		P = self.np_points_
		(N, d) = np.shape(P)
			# Dimension, its 3D, each point has x,y,z values
		d = float(d)
			# add 1 row of all 1's, thats why we use vstack and np.ones
		Q = np.vstack([np.copy(P.T), np.ones(N)]) 
		QT = Q.T

		err = 1.0 + tolerance
		u = (1.0 / N) * np.ones(N)
		
		# Khachiyan Algorithm
		while err > tolerance:
			V = np.dot(Q, np.dot(np.diag(u), QT))

			try:
				D = linalg.inv(V)
			except:
				fix_flag=True
				return (True,fix_flag)
			M = np.diag(np.dot(QT , np.dot(D, Q)))   	
			#M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))
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
		try:
		# the A matrix for the ellipse
			A = linalg.inv(
						   np.dot(P.T, np.dot(np.diag(u), P)) - 
						   np.array([[a * b for b in center] for a in center])
						   ) / d
		except:
			fix_flag=True
			return(True,fix_flag)

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
		#print(self.layers_array)
		bpy.ops.mesh.primitive_uv_sphere_add(size=0.03, location = (0.0,0.0,0.0),layers=self.layers_array)
		bpy.context.object.name = 'fake' + label
		#print(self.layer_indx)                                                                                                                            
		bpy.context.scene.layers[self.layer_index] = True #11 but 12 in count
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
			#print(self.layers_array)
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
			#bpy.data.objects['fake'+label].hide = True
			# we will delete fake objects instead of hiding
			bpy.data.objects['fake'+label].select = True
			bpy.ops.object.delete()

			return False,False
		return False,False


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

	def relationship_lines_off(self):
		for area in bpy.context.screen.areas:
			print(area.type)
			if area.type == 'VIEW_3D':
				print('relationship lines are ',area.spaces[0].show_relationship_lines, 'will be switched off now')
				if area.spaces[0].show_relationship_lines:
					area.spaces[0].show_relationship_lines = False

	#def init(self):
	#	bpy.types.Scene.data_glyc_clusters = []
		
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
		bpy.types.Scene.dict_nobj_sa_vol = {}
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
		if 'Endothelial cell ' in bpy.context.scene.data_names:
			patterns.append('Endothelial cell')
		elif 'endothelial cell' in bpy.context.scene.data_names:
			patterns.append('endothelial cell')
		if 'Pericyte' in bpy.context.scene.data_names:
			patterns.append('Pericyte')
		elif 'pericyte' in bpy.context.scene.data_names:
			patterns.append('pericyte')

		firstTime = True
		if patterns:
			for patt in patterns:
				print("getting ",bpy.context.scene.prop_vertType," for: ", patt)
				#temp_attrib_np, temp_verts_np = getVertices(patt, "Centroids")
				#temp_attrib_np, temp_verts_np = getVertices(patt, "All Vertices") #March 26th, 14:33
				temp_attrib_np, temp_verts_np = getVertices(patt, bpy.context.scene.prop_vertType)
				
				
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
					np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,4],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,0],
					bpy.types.Scene.neur_obj_attrib_np[objctsIndx,1],
					dist
					])
				#create dictionary of:{key: <neuro objects>, value:<surface area & volume >} Join(parent&child)
				if bpy.types.Scene.neur_obj_attrib_np.shape[1]== 4:
					nobj_par_child= " ".join((
						bpy.types.Scene.neur_obj_attrib_np[objctsIndx,0],
						bpy.types.Scene.neur_obj_attrib_np[objctsIndx,1]
						))
					nobj_sa_vol= " ".join((
						bpy.types.Scene.neur_obj_attrib_np[objctsIndx,2],
						bpy.types.Scene.neur_obj_attrib_np[objctsIndx,3]
						))

					bpy.types.Scene.dict_nobj_sa_vol[nobj_par_child] = nobj_sa_vol

				print(
					[np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,0],
					np.array(bpy.types.Scene.data_clusters_centroids)[clusterId,4],
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
	

class OBJECTS_OT_export_clusters_measures(bpy.types.Operator):#tryout
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
			#bpy.types.Scene.dict_nobj_sa_vol[nobj_par_child] = nobj_sa_vol
			print("writing data of [clusters centroids to spine-bouton centroids] distances")
			with open(the_name,'wt') as output:
				writer = csv.writer(output, delimiter='\t')
				writer.writerow(['CId','#granules','ObjectName','ObjectParent','surface area','volume', 'Distance'])
				
				for CId, CCentroid,ObjectName,ObjectParent, Distance in bpy.types.Scene.data_clusters_distances:
					n = ObjectName.rsplit('_',1)
					if bpy.types.Scene.dict_nobj_sa_vol[' '.join((ObjectName,ObjectParent))]:
						surfAreaVolume = bpy.types.Scene.dict_nobj_sa_vol[' '.join((ObjectName,ObjectParent))].rsplit(' ',1)
					else:
						surfAreaVolume=[0,0]
						print('error: no parent-child key found in dict_nobj_sa_vol')
					if not n:
						n[0]='Error splitting name'
					writer.writerow([CId
						, CCentroid
						,n[0]
						,ObjectParent
						,surfAreaVolume[0]
						,surfAreaVolume[1]
						, Distance])
		elif bpy.context.scene.prop_bool_glycogen == True:
			print("writing data of [glycogen to closests neighbors objects- distances]")
			with open(the_name,'wt') as output:
				writer = csv.writer(output, delimiter='\t')
				#writer.writerow(['Neural Object','Parent','No.Associated Glycogens','Glycogen Names','Distance','SurfArea','Volume','Size'])
				writer.writerow(['Neural Object'
					,' --- Parent---'
					, '---SurfArea---'
					,'---Volume---'
					,'----No.Associated Glycogens----'
					,'----Glycogen Names----'
					,'----Size----'
					,'----Distance----',])

				#rowToWrite = []
				for neurObj, glyList in bpy.types.Scene.neuro_glyList_dict.items():
					child_parent = neurObj.rsplit(' ',1)
					child = child_parent[0].rsplit('_',1) #taking out surf*,vol*,solid* strings from names
					if bpy.types.Scene.dict_temp3[neurObj]:
						sa_vol = bpy.types.Scene.dict_temp3[neurObj].rsplit(' ',1)
						surf_area = sa_vol[0]
						volume = sa_vol[1]
					else:
						surf_area = "error in bpy.types.Scene.dict_temp3, neurnal object not found"
					if not child:
						child[0] = 'error splitting object name'
					writeOnce = True
					
					for glyName in glyList:
						gdistance = [dist for gname, dist in bpy.types.Scene.data_glyc_to_neighb_dist if gname == glyName]
						if bpy.types.Scene.dict_temp2[glyName]:
							gly_size = bpy.types.Scene.dict_temp2[glyName]
						else:
							print('error in bpy.types.Scene.dict_temp2, glyname not found in dictionary')
						if gdistance:
							if writeOnce:
								writer.writerow([child[0], child_parent[1], surf_area, volume, len(glyList), glyName, gly_size, gdistance[0]])# gidstance[0]=0.091054856874, gdistance=[0.091054856874324699]
								writeOnce = False
							else:
								writer.writerow([child[0], child_parent[1], surf_area, volume, 0, glyName, gly_size, gdistance[0]])
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
	bl_space_type = "VIEW_3D" 
	bl_region_type =  "UI"

	def draw(self,context):
		layout = self.layout
		scene = context.scene
		obj = context.object
		row1 = layout.row()
		row1.alignment = 'LEFT'
		row1.operator("view3d.activate_addon_button", "Activate/Reset addon")

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
		if bpy.types.Scene.prop_bool_glycogen: #not None or if bpy.types.Scene.prop_bool_clusters
			row1_col1.prop(scene,'prop_bool_glycogen')
			row1_col1.prop(scene,'prop_bool_clusters')

			if bpy.context.scene.prop_bool_clusters == True: 
				row2_layout = layout.row()
				row2_layout.alignment = 'LEFT'
				row2_layout.operator("VIEW3D_OT_show_info", icon = "INFO")

				if bpy.context.scene.prop_dbscan_info:
					infobox = layout.box()
					row2_layout= infobox.row()
					row2_layout.label("Clustering parameters can either be set automatically")
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

				## add type of vertices here --:
				row7_1_step2 = layout.row()
				row7_1_step2.prop(scene,"prop_vertType")
				row7_1_step2.enabled = False
				#
				row8_step2 = layout.row()
				row8_step2.operator("objects.clusters_nearest_neighbours", "Calculate Nearest Neighbour")
				row8_step2.enabled = False

				row9_step2 = layout.row()
				row9_step2.operator("objects.export_clusters_measures", "Export Measurements to File", icon = "DISK_DRIVE")
				row9_step2.enabled = False
				
				if bpy.context.scene.flag_clusters_measure:
					row8_step2.enabled = True
					row7_1_step2.enabled = True
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
					row1_glyc.label("---Output View1---")

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

					# ++++ ALTERNATIVE VIEW +++++ # bouton1 Axon138_E
					row_UIList = layout.row()
					row_UIList.label("---Output View2---")
					layout.template_list("SCENE_UI_List", "", scene, "UIList_glyc_neighb", scene, "UIList_glyc_neighb_indx" )
					if scene.UIList_glyc_neighb_indx >= 0 and len(scene.UIList_glyc_neighb) > 0:
						item = scene.UIList_glyc_neighb[scene.UIList_glyc_neighb_indx]
						UIlistUpdate(item)
						
						col_UIlist = layout.column()
						col_UIlist.prop(item, "li_glyc_neighbours")
						col_UIlist.prop(item, "li_total_granules")
						col_UIlist.enabled = False
						
						row_UIList1 = layout.row()
						#if not bpy.types.Scene.error_flag:
						#row_UIList1.operator("fix.operator","Associated Granules")
						layout.prop(scene,'li_associated_glyco')
						col_UIlist2 = layout.column()
						col_UIlist2.prop(scene,'li_glyc_to_neighb_dist')
						col_UIlist2.enabled = False
					
#--------------------------------------------------------------------------------
#                                   MAIN EXECUTE
#---------------------------------------------------------------------------------
def register():
	
	bpy.utils.register_module(__name__)
	bpy.types.Scene.UIList_glyc_neighb = CollectionProperty(type= PG_List_Entry)#name of property group class
	bpy.types.Scene.UIList_glyc_neighb_indx = IntProperty(name="index for UIList_glyc_neighb", default=0)
	
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
	bpy.types.Scene.prop_min_samples_s2 = None 
	bpy.types.Scene.prop_optimum_eps = None
	bpy.types.Scene.data_glyc_clusters = None
	bpy.types.Scene.prop_nclusters = None
	bpy.types.Scene.prop_silh = None
	#clusters measurements:
	bpy.types.Scene.prop_vertType = bpy.props.EnumProperty(name='Measurements Type:',items=[('Centroids','Centroids',''),('All Vertices','All Vertices','')]) # verts type (cetnroids or all_verts) -constant
	bpy.types.Scene.data_clusters_centroids = None
	bpy.types.Scene.flag_clusters_measure = None
	bpy.types.Scene.flag_clusters_measured = None
	bpy.types.Scene.data_clusters_distances = None
	#glycogens measurements:
	bpy.types.Scene.data_glyc_distances = None
	bpy.types.Scene.data_glyc_distances_occur = None
	bpy.types.Scene.glycogen_attrib_np = np.array([])
	bpy.types.Scene.glycogen_verts_np = np.array([])
	bpy.types.Scene.neur_obj_verts_np = np.array([])
	bpy.types.Scene.neur_obj_attrib_np = np.array([])

	bpy.types.Scene.neuro_gly_glyFreq_dic_sorted = None
	bpy.types.Scene.neuro_glyList_dict = None

	bpy.types.Scene.prop_glyc_neighbours = None
	bpy.types.Scene.prop_associated_glyco = None
	bpy.types.Scene.prop_total_granules = None
	bpy.types.Scene.prop_glyc_to_neighb_dist = None
	
	bpy.types.Scene.data_glyc_neighbours = None
	bpy.types.Scene.data_noOfGlyc = None
	bpy.types.Scene.data_glyc_to_neighb_dist = None

	
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
	if bpy.types.Scene.prop_vertType: #dropdown list(centroids or all vertices)
		del bpy.types.Scene.prop_vertType
	if bpy.types.Scene.data_clusters_centroids:
		del bpy.types.Scene.data_clusters_centroids
	if bpy.types.Scene.flag_clusters_measure:
		del bpy.types.Scene.flag_clusters_measure
	if bpy.types.Scene.flag_clusters_measured:
		del bpy.types.Scene.flag_clusters_measured
	if bpy.types.Scene.data_clusters_distances:
		del bpy.types.Scene.data_clusters_distances
	#glycogen measurements:
	if bpy.types.Scene.data_glyc_distances:
		del bpy.types.Scene.data_glyc_distances 
	if bpy.types.Scene.data_glyc_distances_occur:
		del bpy.types.Scene.data_glyc_distances_occur

	if bpy.types.Scene.glycogen_attrib_np.size>0:
		del bpy.types.Scene.glycogen_attrib_np
	if bpy.types.Scene.glycogen_verts_np.size>0:
		del bpy.types.Scene.glycogen_verts_np
	if bpy.types.Scene.neur_obj_verts_np.size>0:
		del bpy.types.Scene.neur_obj_verts_np
	if bpy.types.Scene.neur_obj_attrib_np.size>0:
		del bpy.types.Scene.neur_obj_attrib_np
	
	if bpy.types.Scene.neuro_gly_glyFreq_dic_sorted:
		del bpy.types.Scene.neuro_gly_glyFreq_dic_sorted
	if bpy.types.Scene.neuro_glyList_dict:
		del bpy.types.Scene.neuro_glyList_dict

	if bpy.types.Scene.prop_glyc_neighbours:
		del bpy.types.Scene.prop_glyc_neighbours
	if bpy.types.Scene.prop_associated_glyco:
		del bpy.types.Scene.prop_associated_glyco
	if bpy.types.Scene.prop_total_granules:
		del bpy.types.Scene.prop_total_granules
	if bpy.types.Scene.prop_glyc_to_neighb_dist:
		del bpy.types.Scene.prop_glyc_to_neighb_dist
	
	if bpy.types.Scene.data_glyc_neighbours:
		del bpy.types.Scene.data_glyc_neighbours
	if bpy.types.Scene.data_noOfGlyc:
		del  bpy.types.Scene.data_noOfGlyc
	if bpy.types.Scene.data_glyc_to_neighb_dist:
		del bpy.types.Scene.data_glyc_to_neighb_dist

	#free memory from UIList gadeget:
	if bpy.types.Scene.UIList_glyc_neighb:
		index = len(bpy.context.scene.UIList_glyc_neighb)
		for i in range(len(bpy.context.scene.UIList_glyc_neighb)+1): #includes 0 index
			list = bpy.context.scene.UIList_glyc_neighb
			list.remove(index)			

			if index >= 0:
				index = index - 1
			else:
				break

	
	
	