class_name RecordingMenu
extends Control
## The recordings found in a directory (recordings/ by default), one card each: a world of Origin or the biphasic
## vessel, with how long it runs and how many frames it holds. A click opens it.

signal chosen(path: String)

var directory: String = ""
var _list: GridContainer


func _init(recordings: String) -> void:
	directory = recordings


func _ready() -> void:
	set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	var background := ColorRect.new()
	background.color = Look.BACKGROUND
	background.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	add_child(background)
	var page := VBoxContainer.new()
	page.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	page.offset_left = 64
	page.offset_top = 48
	page.offset_right = -64
	page.offset_bottom = -48
	page.add_theme_constant_override("separation", 10)
	add_child(page)
	page.add_child(Look.label("protolife — simulations in three dimensions", 34, Look.INK, "title"))
	page.add_child(Look.label("Recordings of simulations, shown as they were computed: nothing is simulated here. "
		+ "A recording in progress is followed live.", 17, Look.MUTED))
	var bar := HBoxContainer.new()
	bar.add_theme_constant_override("separation", 12)
	page.add_child(bar)
	var where := Look.label("Folder: " + directory, 13, Look.MUTED, "mono")
	where.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	bar.add_child(where)
	var refresh := Look.button("Refresh")
	refresh.pressed.connect(fill)
	bar.add_child(refresh)
	var scroll := ScrollContainer.new()
	scroll.size_flags_vertical = Control.SIZE_EXPAND_FILL
	page.add_child(scroll)
	_list = GridContainer.new()
	_list.columns = 3
	_list.add_theme_constant_override("h_separation", 18)
	_list.add_theme_constant_override("v_separation", 18)
	_list.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	scroll.add_child(_list)
	fill()


## Lists the directory again: every subdirectory holding a readable recording.
func fill() -> void:
	for child: Node in _list.get_children():
		child.queue_free()
	var names := DirAccess.get_directories_at(directory)
	names.sort()
	var found := 0
	for name: String in names:
		var recording := Recording.new()
		if recording.open(directory.path_join(name)):
			_list.add_child(_card(recording, name))
			found += 1
	if found == 0:
		_list.add_child(Look.label("No recording here. To make one, from the protolife repository:\n"
			+ "  build/examples/reproduction out=recordings/step-6-random2\n"
			+ "  build/examples/reproduction out=recordings/step-6-control living=0", 15, Look.MUTED, "mono"))


func _card(recording: Recording, name: String) -> Control:
	var card := Button.new()
	card.custom_minimum_size = Vector2(560, 170)
	card.focus_mode = Control.FOCUS_NONE
	var normal := Look.panel(Look.PANEL, 16, 22)
	var hover := normal.duplicate() as StyleBoxFlat
	hover.border_color = Look.PRIMARY
	hover.set_border_width_all(2)
	card.add_theme_stylebox_override("normal", normal)
	card.add_theme_stylebox_override("hover", hover)
	card.add_theme_stylebox_override("pressed", hover)
	card.pressed.connect(func() -> void: chosen.emit(recording.directory))
	var column := VBoxContainer.new()
	column.set_anchors_and_offsets_preset(Control.PRESET_FULL_RECT)
	column.offset_left = 22
	column.offset_top = 18
	column.offset_right = -22
	column.offset_bottom = -18
	column.mouse_filter = Control.MOUSE_FILTER_IGNORE
	column.add_theme_constant_override("separation", 6)
	card.add_child(column)
	var kinds := {"vessel": "VESSEL", "vesicles": "VESICLES", "lineage": "VESICLE LINEAGE", "world": "ORIGIN WORLD", "capsule": "PROTOCELL"}
	var kind := str(kinds.get(recording.kind, recording.kind.to_upper()))
	column.add_child(Look.label(kind + " · " + name, 12, Look.PRIMARY if recording.kind == "world" else Look.LIVING_TEXT, "bold"))
	column.add_child(Look.label(recording.title, 22, Look.INK, "title"))
	var subtitle := Look.label(recording.subtitle, 14, Look.MUTED)
	subtitle.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	column.add_child(subtitle)
	var count := recording.frame_count()
	var last := recording.frame(count - 1)
	var span := Look.time(last.time, recording.time_unit) if last != null else "—"
	column.add_child(Look.label("%d frames · %s simulated" % [count, span], 13, Look.MUTED, "medium"))
	for child: Node in column.get_children():
		(child as Control).mouse_filter = Control.MOUSE_FILTER_IGNORE
	return card
