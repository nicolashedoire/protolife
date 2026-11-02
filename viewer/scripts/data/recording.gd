class_name Recording
extends RefCounted
## A recording of a simulation in the neutral scene format (include/protolife/scene.h): what the manifest says the
## world is and what each number means, and the frames of frames.bin, read as they are appended so that a run can be
## followed live. A frame is the time (float64), the values (float32 each), then the fields (float32, field after
## field over the cells, row-major), little-endian.

const FORMATS := ["protolife-scene-1", "genesis-scene-1"]  ## the same format; Genesis's name too

var directory: String = ""
var kind: String = ""
var title: String = ""
var subtitle: String = ""
var source: String = ""
var time_unit: String = "s"
var n: int = 1
var box: float = 0.0
var frame_bytes: int = 0
var values: Array[Dictionary] = []
var fields: Array[Dictionary] = []
var constants: Dictionary = {}
var constant_labels: Dictionary = {}

var _value_index: Dictionary = {}
var _field_index: Dictionary = {}
var _file: FileAccess


## Reads the manifest of a recording's directory; false if it is not one.
func open(path: String) -> bool:
	directory = path
	var text := FileAccess.get_file_as_string(path.path_join("scene.json"))
	if text.is_empty():
		return false
	var parsed: Variant = JSON.parse_string(text)
	if not parsed is Dictionary:
		return false
	var manifest: Dictionary = parsed
	if not FORMATS.has(str(manifest.get("format", ""))):
		return false
	kind = str(manifest.get("kind", ""))
	title = str(manifest.get("title", ""))
	subtitle = str(manifest.get("subtitle", ""))
	source = str(manifest.get("source", ""))
	time_unit = str(manifest.get("time_unit", "s"))
	frame_bytes = int(manifest.get("frame_bytes", 0))
	var grid: Dictionary = manifest.get("grid", {})
	n = int(grid.get("n", 1))
	box = float(grid.get("box_m", 0.0))
	values = _quantities(manifest, "values")
	fields = _quantities(manifest, "fields")
	for i in values.size():
		_value_index[str(values[i]["key"])] = i
	for i in fields.size():
		_field_index[str(fields[i]["key"])] = i
	for constant: Dictionary in _quantities(manifest, "constants"):
		constants[str(constant["key"])] = float(constant.get("value", 0.0))
		constant_labels[str(constant["key"])] = str(constant.get("label", ""))
	var expected := 8 + 4 * (values.size() + fields.size() * n * n)
	return frame_bytes == expected


func _quantities(manifest: Dictionary, name: String) -> Array[Dictionary]:
	var list: Array[Dictionary] = []
	var raw: Variant = manifest.get(name, [])
	if raw is Array:
		for item: Variant in raw:
			if item is Dictionary:
				list.append(item)
	return list


## The whole frames the file holds now; it grows while the run goes on.
func frame_count() -> int:
	if _file == null:
		_file = FileAccess.open(directory.path_join("frames.bin"), FileAccess.READ)
	if _file == null or frame_bytes <= 0:
		return 0
	return int(_file.get_length() / frame_bytes)


## One frame, or null when it is not there yet.
func frame(index: int) -> RecordingFrame:
	if index < 0 or index >= frame_count():
		return null
	_file.seek(index * frame_bytes)
	var result := RecordingFrame.new()
	result.index = index
	result.time = _file.get_double()
	if not values.is_empty():
		result.values = _file.get_buffer(4 * values.size()).to_float32_array()
	if not fields.is_empty():
		result.fields = _file.get_buffer(4 * fields.size() * n * n).to_float32_array()
	return result


func value_index(key: String) -> int:
	return int(_value_index.get(key, -1))


func field_index(key: String) -> int:
	return int(_field_index.get(key, -1))


## A value of a frame by its key; NAN when the recording has no such value.
func value(of: RecordingFrame, key: String) -> float:
	var i := value_index(key)
	return of.values[i] if i >= 0 and i < of.values.size() else NAN


func constant(key: String, fallback: float = NAN) -> float:
	return float(constants.get(key, fallback))
