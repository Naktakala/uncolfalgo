print("Executing input file")
--############################################### Mesh
chiMeshHandlerCreate()

Lx = 5.0
Ly = 5.0
Nx = 40
Ny = 40

dx = Lx/Nx
dy = Ly/Ny
nodesx={}
nodesy={}
for i=0,Nx do
    nodesx[i+1] = i*dx
end
for j=0,Ny do
    nodesy[j+1] = j*dy
end

chiMeshCreate2DOrthoMesh(nodesx,nodesy)
chiVolumeMesherExecute();

--########################################## Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--vol1 = chiLogicalVolumeCreate(RPP,2.0,2.5,2.0,2.5,-1000,1000)
vol1 = chiLogicalVolumeCreate(RPP,0.0,0.5*10/Nx,0.0,0.5*10/Ny,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--########################################## Export mesh
print("Exporting mesh")
chiRegionExportMeshToVTK(0,"YMesh")