# MDHTTracking
Multi-dimensional Hough Transform tracking

### Getting started:

To complie the library:
```
cmake -S MDHTTracking -B build
cmake --build build
```
Steering file configuration:
```
HTATree = MarlinProcessorWrapper("HTATree")
HTATree.OutputLevel = INFO
HTATree.ProcessorType = "HTATrainingTree"
HTATree.Parameters = {
    "MCParticleCollection": [ MCP ],
    "TrackerHitCollections": [ VXDB_d, VXDE_d, ITDB_d, ITDE_d, OTDB_d, OTDE_d ],
    "TrackerSimHitCollections": [ VXDB_s, VXDE_s, ITDB_s, ITDE_s, OTDB_s, OTDE_s ],
    "TrackerHitRelationCollections": [ VXDB_r, VXDE_r, ITDB_r, ITDE_r, OTDB_r, OTDE_r ],
    "ParticleTypes": [ "13" ],
    "SaveOnlyPartHits": [ "true" ]
}
```
