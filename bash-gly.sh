#!/bin/bash
if [[ "$OSTYPE" == "linux-gnu" ]]; then
	if [ -d "/usr/local/lib/python3.4:/usr/local/lib/python3.4/site-packages" ] && [ -d "/usr/local/lib/python3.4/lib-dynload" ]; then
		export PYTHONPATH=/usr/local/lib/python3.4:/usr/local/lib/python3.4/site-packages:/usr/local/lib/python3.4/lib-dynload:$PYTHONPATH
	elif [ -d "$HOME/anaconda/lib/python3.4/site-packages" ]; then
		export PYTHONPATH=$HOME/anaconda/lib/python3.4/site-packages:$PYTHONPATH
	fi


elif [[ "$OSTYPE" == "darwin"* ]]; then
	if [[ $(sysctl -n hw.ncpu) > "1" ]]; then
		echo "$(sysctl -n hw.ncpu)"
		#let "TOT=$(sysctl -n hw.ncpu)/2"
		let "TOT=$(sysctl -n hw.ncpu)"
		export OMP_NUM_THREADS=$TOT
	fi

	if [ -d "$HOME/anaconda/lib/python3.4/site-packages" ]; then
		export PYTHONPATH=$HOME/anaconda/lib/python3.4/site-packages:$PYTHONPATH
		
	elif [ -d "/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages" ] && [ -d "/usr/local/bin/python3.4" ]; then
		export PATH=/usr/local/bin/python3.4:$PATH
		export PYTHONPATH=/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages:$PYTHONPATH
	fi
	
fi

if [ -d "/Applications/blender.app/Contents/MacOS" ]; then
	echo "blender is on"
	/Applications/blender.app/Contents/MacOS/blender
else
	echo "Blender not found in Applications. please modify path"
fi


