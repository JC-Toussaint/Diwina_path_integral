# FeeLLGood micromagnetic imaging toolbox (fmit)

fmit-GUI :
python executable that launch global application with a GUI, to build some settings, save those
settings in yaml with some 3D visualization for simulating phase images and XMCD images
Compute button launch a C++ executable (pathIntegral)

phase_contours.py : iso values of the holographic phase

fresnel_from_png_phase.py : post-processing : computes Fresnel intensities (to check)

fg-yml.py : 'tree like' text display of a yaml file 

msh-explorer.py : informations on a mesh file: number of triangles, tetrahedrons, ...

fmit-png-extractor : tool to get some metadatas from an png file from pathIntegral

fmit-png-viewer : tool to display pngs output from pathIntegral

