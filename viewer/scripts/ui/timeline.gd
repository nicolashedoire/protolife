class_name Timeline
extends PanelContainer
## The time bar under every view: play or pause, the frames of the recording to scrub, the playing speed in frames
## per second, and "en direct", which follows the last frame while a run is still writing.

signal frame_changed(index: int)

const SPEEDS := [2.0, 8.0, 30.0, 90.0]
const SPEED_NAMES := ["2 frames/s", "8 frames/s", "30 frames/s", "90 frames/s"]

var recording: Recording
var current: int = 0
var playing: bool = true
var live: bool = false
var _speed: float = 8.0
var _position: float = 0.0
var _count: int = 0
var _poll: float = 0.0
var _play: Button
var _slider: HSlider
var _time: Label
var _live: Button
var _speeds: Array[Button] = []


func _init(source: Recording) -> void:
	recording = source


func _ready() -> void:
	add_theme_stylebox_override("panel", Look.panel(Look.PANEL, 0, 14))
	set_anchors_and_offsets_preset(Control.PRESET_BOTTOM_WIDE)
	offset_top = -72
	var row := HBoxContainer.new()
	row.add_theme_constant_override("separation", 14)
	add_child(row)
	_play = Look.button("Pause", true)
	_play.custom_minimum_size = Vector2(130, 0)
	_play.pressed.connect(func() -> void: set_playing(not playing))
	row.add_child(_play)
	_slider = HSlider.new()
	_slider.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_slider.size_flags_vertical = Control.SIZE_SHRINK_CENTER
	_slider.step = 1.0
	_slider.value_changed.connect(_on_slider)
	row.add_child(_slider)
	_time = Look.label("", 15, Look.INK, "medium")
	_time.custom_minimum_size = Vector2(250, 0)
	row.add_child(_time)
	for i in SPEEDS.size():
		var speed := Look.button(str(SPEED_NAMES[i]), float(SPEEDS[i]) == _speed)
		speed.pressed.connect(_on_speed.bind(i))
		_speeds.append(speed)
		row.add_child(speed)
	_live = Look.button("Live", live)
	_live.pressed.connect(func() -> void: set_live(not live))
	row.add_child(_live)
	_count = recording.frame_count()
	_refresh_slider()
	_show(0)


func set_playing(value: bool) -> void:
	playing = value
	_play.text = "Pause" if playing else "Play"


func set_live(value: bool) -> void:
	live = value
	_restyle(_live, live)
	if live:
		_show(_count - 1)


## Shows one frame (clamped), and tells the view.
func show_frame(index: int) -> void:
	_show(index)


func _process(delta: float) -> void:
	_poll += delta
	if _poll >= 0.5:  # a run still writing grows its file
		_poll = 0.0
		var count := recording.frame_count()
		if count != _count:
			_count = count
			_refresh_slider()
			if live:
				_show(_count - 1)
	if playing and not live and _count > 1:
		_position += _speed * delta
		var next := int(_position)
		if next >= _count:
			next = _count - 1
			set_playing(false)
		if next != current:
			_show(next)


func _show(index: int) -> void:
	current = clampi(index, 0, maxi(_count - 1, 0))
	_position = float(current)
	_slider.set_value_no_signal(current)
	var frame := recording.frame(current)
	var total := recording.frame(_count - 1)
	_time.text = "%s / %s" % [Look.time(frame.time if frame else NAN, recording.time_unit), Look.time(total.time if total else NAN, recording.time_unit)]
	frame_changed.emit(current)


func _refresh_slider() -> void:
	_slider.max_value = maxf(_count - 1, 0)


func _on_slider(value: float) -> void:
	set_live(false)
	_show(int(value))


func _on_speed(which: int) -> void:
	_speed = float(SPEEDS[which])
	for i in _speeds.size():
		_restyle(_speeds[i], i == which)


func _restyle(button: Button, active: bool) -> void:
	var fresh := Look.button(button.text, active)
	for style: String in ["normal", "hover", "pressed"]:
		button.add_theme_stylebox_override(style, fresh.get_theme_stylebox(style))
	for colour: String in ["font_color", "font_hover_color", "font_pressed_color"]:
		button.add_theme_color_override(colour, fresh.get_theme_color(colour))
	fresh.free()
