PARAVIEW=paraview
if ! [ -x "$(command -v paraview)" ]; then
	PARAVIEW=/Applications/ParaView-5.4.1.app/Contents/MacOS/paraview
fi
