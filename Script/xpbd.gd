extends Node3D

var tet_mesh: TetGenMesh
var tet_array_constructor: TetArrayConstructor
var mesh_instance: MeshInstance3D
var surface_mesh: ArrayMesh

var pos: PackedVector3Array
var prev_pos: PackedVector3Array
var inv_mass: PackedFloat32Array
var velocity: PackedVector3Array

var edge_Ids: PackedInt32Array
var tet_Ids: PackedInt32Array

var edge_compliance: float = 0.0001
var volume_compliance: float = 0.000001

var edge_lengths: PackedFloat32Array
var tet_volumes: PackedFloat32Array

# [e1, e2, e3 | e2, e3 | e3, e1] Compact Edge Array
# [2, 1, 1] Edge Stride
# [0, 3, 5] Edge Vertex Indices
func solveEdges(compliance: float = 0.0001, dt: float = 0.01):
	var alpha: float = compliance / dt / dt
	for i in range(0, edge_Ids.size(), 2):
		var id1: int = edge_Ids[i]
		var id2: int = edge_Ids[i+1]
		var v1: Vector3 = pos[id1]
		var v2: Vector3 = pos[id2]
		var w1: float = inv_mass[id1]
		var w2: float = inv_mass[id2]
		var w: float = w1 + w2
		if w == 0.0:
			continue
		
		var grad: Vector3 = v1 - v2
		var len: float = grad.length()
		
		if len == 0.0:
			continue
		
		grad = grad.normalized()
		var rest_len: float = edge_lengths[i]
		var C: float = len - rest_len
		var s: float = -C / (w + alpha)
		pos[id1] += grad * s * w1
		pos[id2] += grad * -s * w2

# [e1, e2, e3, e4 | e2, e6, e10, e1, e11, e20, e12] Compact Tetrahedron Array
# [3, 6] # Tetrahedron Stride | Always multiple of 3
# [0, 4] # Tetrahedron Vertex Indices
func solveVolumes(compliance: float, dt: float):
	var alpha: float = compliance / dt / dt
	var volIdOrder = [[1,2,3], [0,2,3], [0,3,1], [0,1,2]]  # Same as his volIdOrder

	for i in range(0, tet_Ids.size(), 4):
		var w = 0.0
		var grads = []  # Store gradients for each vertex
		grads.resize(4)
		
		# Compute gradients for all 4 vertices using his ordering
		for j in range(4):
			var id0 = tet_Ids[i + volIdOrder[j][0]]
			var id1 = tet_Ids[i + volIdOrder[j][1]]
			var id2 = tet_Ids[i + volIdOrder[j][2]]
			
			# Get positions
			var p0 = pos[id0]
			var p1 = pos[id1]
			var p2 = pos[id2]
			
			# Compute differences
			var v1 = p1 - p0
			var v2 = p2 - p0
			
			# Cross product and scale by 1/6
			var grad = v1.cross(v2) * (1.0/6.0)
			grads[j] = grad
			
			# Accumulate weighted sum
			w += inv_mass[tet_Ids[i + j]] * grad.length_squared()
		
		if w == 0.0:
			continue
		
		# Compute volume and constraint
		var vol = getTetVolume(i)  # Make sure this uses the tetrahedron at index i
		var rest_vol = tet_volumes[i]
		var C = vol - rest_vol
		var s = -C / (w + alpha)
		
		# Update all four vertices
		for j in range(4):
			var id = tet_Ids[i + j]
			pos[id] += grads[j] * s * inv_mass[id]

func getTetVolume(base_index: int) -> float:
	var base: int = base_index
	var id1: int = tet_Ids[base]
	var id2: int = tet_Ids[base + 1]
	var id3: int = tet_Ids[base + 2]
	var id4: int = tet_Ids[base + 3]
	
	var vec_21: Vector3 = pos[id2] - pos[id1]
	var vec_31: Vector3 = pos[id3] - pos[id1]
	var vec_41: Vector3 = pos[id4] - pos[id1]
	
	return vec_41.dot(vec_21.cross(vec_31)) / 6.0

func computeEdgeRestLengths():
	for i in range(0, edge_Ids.size(), 2):
		var id1: int = edge_Ids[i]
		var id2: int = edge_Ids[i + 1]
		var vert1: Vector3 = pos[id1]
		var vert2: Vector3 = pos[id2]
		var len: float = (vert1 - vert2).length()
		edge_lengths[i] = len

func computeTetRestVolumes():
	for i in range(0, tet_Ids.size(), 4):
		var id1: int = tet_Ids[i]
		var id2: int = tet_Ids[i + 1]
		var id3: int = tet_Ids[i + 2]
		var id4: int = tet_Ids[i + 3]
			
		var vol: float = getTetVolume(i)
		tet_volumes[i] = vol

func preSolve(dt: float, force: Vector3):
	for i in range(pos.size()):
		if (inv_mass[i] == 0.0):
			continue
		velocity[i] += force * dt
		prev_pos[i] = pos[i]
		pos[i] += velocity[i] * dt
		
		if pos[i].y > 3:
			pos[i] = prev_pos[i]
			pos[i].y = 3

func solve(dt: float):
	solveEdges(edge_compliance, dt)
	solveVolumes(volume_compliance, dt)
	#print("Position: ", pos[0])

func postSolve(dt):
	for i in range(pos.size()):
		if inv_mass[i] == 0.0:
			continue
		velocity[i] = (pos[i] - prev_pos[i]) * (1.0 / dt)
	# update mesh
	surface_mesh.clear_surfaces()
	surface_mesh = tet_mesh.update_mesh(pos)
	mesh_instance.mesh = surface_mesh

func _ready():
	tet_mesh = TetGenMesh.new()
	tet_mesh.load_from_base_name("res://Mesh/LowPoly/Suzanne")
	tet_array_constructor = TetArrayConstructor.new()
	
	# Create surface mesh
	surface_mesh = tet_mesh.create_surface_mesh()
	if surface_mesh:
		mesh_instance = MeshInstance3D.new()
		mesh_instance.mesh = surface_mesh
		
		# Center the mesh
		var bounds = surface_mesh.get_aabb()
		var center = bounds.get_center()
		mesh_instance.position = Vector3(0, 2, 0)
		
		# Add to scene
		add_child(mesh_instance)
		
		print("Mesh added to scene. Bounds: ", bounds)
		print("Mesh center: ", center)
		
		pos = tet_mesh.vertices
		prev_pos.resize(pos.size())
		inv_mass.resize(pos.size())
		velocity.resize(pos.size())
		prev_pos.fill(Vector3(0, 0, 0))
		inv_mass.fill(0.5)
		velocity.fill(Vector3(0, 0, 0))
		
		edge_Ids = tet_array_constructor.construct_edge(tet_mesh.edges)
		tet_Ids = tet_array_constructor.construct_tet(tet_mesh.tetrahedra)
		edge_lengths.resize(edge_Ids.size())
		tet_volumes.resize(tet_Ids.size())
		tet_volumes.fill(0.0)
		edge_lengths.fill(0.0)
		computeEdgeRestLengths()
		computeTetRestVolumes()
		
		print("Edge Ids: ", edge_Ids.size())
		print("Tet Ids: ", tet_Ids.size())
		print("Edge Rest Lengths: ", edge_lengths.size())
		print("Tet Rest Volumes: ", tet_volumes.size())
		print("Pos: ", pos.size())

func _process(delta: float) -> void:
	var force: Vector3 = Vector3(0, 5, 0)
	var dt: float = 0.01
	var sdt: float = dt / 10
	preSolve(dt, force)
	solve(dt)
	postSolve(dt)
