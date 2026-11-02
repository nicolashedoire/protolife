class_name OrbitCamera
extends Camera3D
## A camera that turns around a point: drag with the left button to turn, the wheel to come closer or move away.
## It turns slowly by itself until the viewer takes it in hand.

var target := Vector3.ZERO
var distance: float = 4.0
var yaw: float = 0.6
var pitch: float = 0.35
var near_limit: float = 1.5
var far_limit: float = 30.0
var drifting: bool = true
var _dragging: bool = false


func _ready() -> void:
	_place()


func _process(delta: float) -> void:
	if drifting:
		yaw += 0.05 * delta
		_place()


func _unhandled_input(event: InputEvent) -> void:
	if event is InputEventMouseButton:
		var button := event as InputEventMouseButton
		if button.button_index == MOUSE_BUTTON_LEFT:
			_dragging = button.pressed
			drifting = false
		elif button.button_index == MOUSE_BUTTON_WHEEL_UP and button.pressed:
			distance = maxf(near_limit, distance * 0.9)
			_place()
		elif button.button_index == MOUSE_BUTTON_WHEEL_DOWN and button.pressed:
			distance = minf(far_limit, distance * 1.1)
			_place()
	elif event is InputEventMouseMotion and _dragging:
		var motion := event as InputEventMouseMotion
		yaw -= motion.relative.x * 0.006
		pitch = clampf(pitch + motion.relative.y * 0.006, -0.2, 1.45)
		_place()


func _place() -> void:
	var offset := Vector3(cos(pitch) * sin(yaw), sin(pitch), cos(pitch) * cos(yaw)) * distance
	position = target + offset
	look_at(target, Vector3.UP)
