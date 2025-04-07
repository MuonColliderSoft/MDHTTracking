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
    "TrackerHitInputCollections": [ VXDB_d, VXDE_d, ITDB_d, ITDE_d, OTDB_d, OTDE_d ],
    "ParticleTypes": [ "13" ]                                                                                                           
}
```
