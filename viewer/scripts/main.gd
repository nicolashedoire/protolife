extends Node
## The viewer's entry: the menu of recordings, then the view of the one chosen, with the time bar and a way back.
## Options after "--" on the command line: --recordings=DIR (recordings/ of the protolife repository by default),
## --open=NAME (opens that recording at once), --frame=N (shows that frame, paused), --field=KEY (a world's field),
## --shot=FILE.png (saves the screen after a moment and quits: for checks and documents).

var _recordings: String = ""
var _options: Dictionary = {}
var _view: Node
var _timeline: Timeline
var _menu: RecordingMenu
var _shot_in: int = -1


func _ready() -> void:
	for argument: String in OS.get_cmdline_user_args():
		var parts := argument.trim_prefix("--").split("=", true, 1)
		_options[parts[0]] = parts[1] if parts.size() > 1 else ""
	var default := ProjectSettings.globalize_path("res://").path_join("../recordings").simplify_path()
	_recordings = str(_options.get("recordings", default))
	if _options.has("open"):
		_open(_recordings.path_join(str(_options["open"])))
	else:
		_show_menu()
	if _options.has("shot"):
		_shot_in = 45


func _process(_delta: float) -> void:
	if _shot_in < 0:
		return
	_shot_in -= 1
	if _shot_in == 0:
		var image := get_viewport().get_texture().get_image()
		image.save_png(str(_options["shot"]))
		get_tree().quit()


func _show_menu() -> void:
	_clear()
	_menu = RecordingMenu.new(_recordings)
	_menu.chosen.connect(_open)
	var layer := CanvasLayer.new()
	layer.add_child(_menu)
	add_child(layer)
	_view = layer


func _open(path: String) -> void:
	var recording := Recording.new()
	if not recording.open(path):
		push_error("not a protolife recording: " + path)
		_show_menu()
		return
	_clear()
	var view: Node
	if recording.kind == "vessel":
		view = VesselView.new(recording)
	elif recording.kind == "vesicles":
		view = VesiclesView.new(recording)
	elif recording.kind == "lineage":
		view = LineageView.new(recording)
	elif recording.kind == "capsule":
		view = CapsuleView.new(recording)
	else:
		view = WorldView.new(recording, str(_options.get("field", "")))
	add_child(view)
	_view = view
	var layer := CanvasLayer.new()
	layer.layer = 10
	view.add_child(layer)
	var back := Look.button("← Recordings")
	back.position = Vector2(40, 24)
	back.pressed.connect(_show_menu)
	layer.add_child(back)
	_timeline = Timeline.new(recording)
	_timeline.frame_changed.connect(func(index: int) -> void: view.call("show_frame", index))
	layer.add_child(_timeline)
	if _options.has("frame"):
		_timeline.set_playing(false)
		_timeline.call_deferred("show_frame", int(_options["frame"]))
		_options.erase("frame")


func _clear() -> void:
	if _view != null:
		_view.queue_free()
	_view = null
