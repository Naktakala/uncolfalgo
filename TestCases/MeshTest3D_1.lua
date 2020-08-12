print("Executing input file")
--############################################### Mesh
chiMeshHandlerCreate()

if (L == nil) then L = 5.0; end
if (N == nil) then N = 40; end

ds = L/N

nodes={}
for i=0,N do
    nodes[i+1] = i*ds
end

chiMeshCreate3DOrthoMesh(nodes,nodes,nodes)
chiVolumeMesherExecute();

--########################################## Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,0.0,0.5*10/N,
                                  0.0,0.5*10/N,
                                  0.0,0.5*10/N)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--########################################## Export mesh
print("Exporting mesh")
chiRegionExportMeshToVTK(0,"YMesh")