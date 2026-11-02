class_name RecordingFrame
extends RefCounted
## One frame of a recording: its index, its time (s), its values and its fields (field after field over the cells).

var index: int = 0
var time: float = 0.0
var values := PackedFloat32Array()
var fields := PackedFloat32Array()


## One field over the cells, by its index in the manifest.
func field(which: int, cells: int) -> PackedFloat32Array:
	return fields.slice(which * cells, (which + 1) * cells)
