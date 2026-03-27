extends Node3D

var tet_mesh: TetGenMesh
var tet_array_constructor: TetArrayConstructor

var pos: PackedVector3Array
var prev_pos: PackedVector3Array
var inv_mass: PackedFloat32Array
var velocity: PackedVector3Array

var edge_Ids: PackedInt32Array
var tet_Ids: PackedInt32Array

var edge_stride: PackedInt32Array
var edge_vertex_indices: PackedInt32Array

var tet_stride: PackedInt32Array
var tet_vertex_indices: PackedInt32Array

var edge_compliance: float = 0.0
var volume_compliance: float = 0.0

var edge_lengths: PackedFloat32Array
var tet_volumes: PackedFloat32Array

var tet_task_id: int
var edge_task_id: int

# [e1, e2, e3 | e2, e3 | e3, e1] Compact Edge Array
# [2, 1, 1] Edge Stride
# [0, 3, 5] Edge Vertex Indices
func solveEdges(i: int, compliance: float = 0.0, dt: float = 0.01):
	var alpha: float = compliance / dt / dt;
	var id1_index: int = edge_vertex_indices[i]
	var id1: int = edge_Ids[id1_index]
	var stride = edge_stride[i]
	for stri in range(1, stride + 1):
		var id2_index: int = id1_index + stri
		var id2: int = edge_Ids[id2_index]
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
		var rest_len: float = edge_lengths[id2_index]
		var C: float = len - rest_len
		var s: float = -C / (w + alpha)
		pos[id1] += grad * s * w1

# [e1, e2, e3, e4 | e2, e6, e10, e1, e11, e20, e12] Compact Tetrahedron Array
# [3, 6] # Tetrahedron Stride | Always multiple of 3
# [0, 4] # Tetrahedron Vertex Indices
func solveVolumes(i: int, compliance: float = 0.0, dt: float = 0.01):
	var alpha: float = compliance / dt / dt
	var id1_index: int = tet_vertex_indices[i]
	var id1: int = tet_Ids[id1_index]
	var stride: float = tet_stride[i]
	
	var num_iter: int = stride / 3
	var tet_neighbours: int = 1
	for iter in range(num_iter):
		var w = 0.0
		var id2: int = tet_Ids[id1_index + tet_neighbours]
		var id3: int = tet_Ids[id1_index + tet_neighbours + 1]
		var id4: int = tet_Ids[id1_index + tet_neighbours + 2]
		# Precompute position differences once
		var p1 = pos[id1]
		var p2 = pos[id2]
		var p3 = pos[id3]
		var p4 = pos[id4]
		
		var v21 = p2 - p1
		var v31 = p3 - p1
		var v41 = p4 - p1
		var v32 = p3 - p2
		var v42 = p4 - p2
		
		# Gradients (each is 1/6 * cross product of two edge vectors)
		var grad1 = v32.cross(v42) * (1.0/6.0)
		var grad2 = v31.cross(v41) * (1.0/6.0)
		var grad3 = v21.cross(v41) * (1.0/6.0)  # Or compute directly
		var grad4 = v21.cross(v31) * (1.0/6.0)
		
		w += (inv_mass[id1] * grad1.length_squared()) + (inv_mass[id2] * grad2.length_squared()) + (inv_mass[id3] * grad3.length_squared()) + (inv_mass[id4] * grad4.length_squared())
		
		if w == 0.0:
			continue
		
		var vol: float = getTetVolume(id1, tet_neighbours)
		var rest_vol: float = tet_volumes[id1_index + tet_neighbours]
		var C: float = vol - rest_vol
		var s: float = -C / (w + alpha)
		pos[id1] += grad1 * s * inv_mass[id1]
		
		tet_neighbours += 3

func getTetVolume(base_index: int, stride: int) -> float:
	var base: int = base_index
	var id1: int = tet_Ids[base]
	var id2: int = tet_Ids[base + stride]
	var id3: int = tet_Ids[base + stride + 1]
	var id4: int = tet_Ids[base + stride + 2]
	
	var vec_21: Vector3 = pos[id2] - pos[id1]
	var vec_31: Vector3 = pos[id3] - pos[id1]
	var vec_41: Vector3 = pos[id4] - pos[id1]
	
	return vec_41.dot(vec_21.cross(vec_31)) / 6.0

func computeEdgeRestLengths():
	for i in range(edge_vertex_indices.size()):
		var id1_index: int = edge_vertex_indices[i]
		var id1: int = edge_Ids[id1_index]
		var s: int = edge_stride[i]
		edge_lengths[id1_index] = 0.0
		for str in range(1, s + 1):
			var id2: int = edge_Ids[id1_index + str]
			var id2_index: int = id1_index + str
			var vert1: Vector3 = pos[id1]
			var vert2: Vector3 = pos[id2]
			var len: float = (vert1 - vert2).length()
			edge_lengths[id2_index] = len

func computeTetRestVolumes():
	for i in range(tet_vertex_indices.size()):
		var id1_index: int = tet_vertex_indices[i]
		var id1: int = tet_Ids[id1_index]
		var stride: float = tet_stride[i]
		tet_volumes[id1_index] = 0.0
		
		var num_iter: int = stride / 3
		var tet_neighbours: int = 1
		for iter in range(num_iter):
			var id2_index: int = id1_index + tet_neighbours
			var id3_index: int = id1_index + tet_neighbours + 1
			var id4_index: int = id1_index + tet_neighbours + 2
			
			var vol: float = getTetVolume(id1_index, tet_neighbours)
			if vol <= 0.0:
				var temp: int = tet_Ids[id3_index]
				tet_Ids[id3_index] = tet_Ids[id2_index]
				tet_Ids[id2_index] = temp
				vol = getTetVolume(id1_index, tet_neighbours)
			
			tet_volumes[id2_index] = vol
			tet_volumes[id3_index] = vol
			tet_volumes[id4_index] = vol
			
			tet_neighbours += 3

func preSolve(dt: float, force: Vector3):
	for i in range(pos.size()):
		if (inv_mass[i] == 0.0):
			continue
		velocity[i] += force * dt
		prev_pos[i] = pos[i]
		pos[i] += velocity[i] * dt
		
		if pos[i].y < 0.0:
			pos[i] = prev_pos[i]
			pos[i].y = 0.0;

func solve(dt: float):
	tet_task_id = WorkerThreadPool.add_group_task(solveVolumes, tet_vertex_indices.size(), -1, true)
	edge_task_id = WorkerThreadPool.add_group_task(solveEdges, edge_vertex_indices.size(), -1, true)
	WorkerThreadPool.wait_for_group_task_completion(tet_task_id)
	WorkerThreadPool.wait_for_group_task_completion(edge_task_id)

func postSolve(dt):
	for i in range(pos.size()):
		if inv_mass[i] == 0.0:
			continue
		velocity[i] = (pos[i] - prev_pos[i]) * (1.0 / dt)
	# update mesh

func _ready():
	tet_mesh = TetGenMesh.new()
	tet_mesh.load_from_base_name("res://Mesh/LowPoly/Suzanne")
	tet_array_constructor = TetArrayConstructor.new()
	
	# Create surface mesh
	var surface_mesh = tet_mesh.create_surface_mesh()
	if surface_mesh:
		var mesh_instance = MeshInstance3D.new()
		mesh_instance.mesh = surface_mesh
		
		# Center the mesh
		var bounds = surface_mesh.get_aabb()
		var center = bounds.get_center()
		mesh_instance.position = -center
		
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
		
		edge_Ids = tet_array_constructor.construct_compact_edge(tet_mesh.edges, edge_stride, edge_vertex_indices)
		tet_Ids = tet_array_constructor.construct_compact_tet(tet_mesh.tetrahedra, tet_stride, tet_vertex_indices)
		edge_lengths.resize(edge_Ids.size())
		tet_volumes.resize(tet_Ids.size())
		tet_volumes.fill(0.0)
		computeEdgeRestLengths()
		computeTetRestVolumes()
		
		print("Edge Ids: ", edge_Ids.size())
		print("Tet Ids: ", tet_Ids.size())
		print("Edge Stride: ", edge_stride.size())
		print("Edge Vertex Indices: ", edge_vertex_indices.size())
		print("Edge Rest Lengths: ", edge_lengths.size())
		print("Tet Stride: ", tet_stride.size())
		print("Tet Vertex Indices: ", tet_vertex_indices.size())
		print("Tet Rest Volumes: ", tet_volumes.slice(0, 20))
		print("Pos: ", pos.size())

func _process(delta: float) -> void:
	var force: Vector3 = Vector3(0, 1, 0)
	preSolve(delta, force)
	solve(delta)
	postSolve(delta)
