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
	solveEdges(edge_compliance, dt)
	solveVolumes(volume_compliance, dt)

# [e1, e2, e3 | e2, e3 | e3, e1] Compact Edge Array
# [2, 1, 1] Edge Stride
# [0, 3, 5] Edge Vertex Indices
func solveEdges(compliance: float, dt: float):
	var alpha: float = compliance / dt / dt;
	for i in range(edge_vertex_indices.size()):
		var id1: int = edge_Ids[edge_vertex_indices[i]]
		var stride = edge_stride[i]
		for stri in range(1, stride + 1):
			var id2: int = edge_Ids[id1 + stri]
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
			var rest_len: float = edge_lengths[id2]
			var C: float = len - rest_len
			var s: float = -C / (w + alpha)
			pos[id1] += grad * s * w1

# [e1, e2, e3, e4 | e2, e6, e10, e1, e11, e20, e12] Compact Tetrahedron Array
# [3, 6] # Tetrahedron Stride | Always multiple of 3
# [0, 4] # Tetrahedron Vertex Indices
func solveVolumes(compliance: float, dt: float):
	var alpha: float = compliance / dt / dt
	for i in range(tet_vertex_indices.size()):
		var w = 0.0
		var id1: int = tet_Ids[tet_vertex_indices[i]]
		var stride: float = tet_stride[i]
		
		var num_iter: int = stride / 3
		var tet_neighbours: int = 1
		for iter in range(num_iter):
			var id2: int = tet_Ids[id1 + tet_neighbours]
			var id3: int = tet_Ids[id2 + tet_neighbours + 1]
			var id4: int = tet_Ids[id3 + tet_neighbours + 2]
			
			# for id1
			var vec_32: Vector3 = pos[id3] - pos[id2]
			var vec_42: Vector3 = pos[id4] - pos[id2]
			var grad1: Vector3 = vec_32.cross(vec_42)
			grad1 /= (1.0 / 6.0)
			
			# ... id2
			var vec_31: Vector3 = pos[id3] - pos[id1]
			var vec_41: Vector3 = pos[id4] - pos[id1]
			var grad2: Vector3 = vec_31.cross(vec_41)
			grad2 /= (1.0 / 6.0)
			
			# ... id3
			var vec_21: Vector3 = pos[id2] - pos[id1]
			vec_41 = pos[id4] - pos[id1]
			var grad3: Vector3 = vec_31.cross(vec_41)
			grad3 /= (1.0 / 6.0)
			
			# ... id4
			vec_21 = pos[id2] - pos[id1]
			vec_31 = pos[id3] - pos[id1]
			var grad4: Vector3 = vec_31.cross(vec_41)
			grad4 /= (1.0 / 6.0)
			
			w += (inv_mass[id1] * grad1.length_squared()) + (inv_mass[id2] * grad2.length_squared()) + (inv_mass[id3] * grad3.length_squared()) + (inv_mass[id4] * grad4.length_squared())
			
			if w == 0.0:
				continue
			
			var vol: float = getTetVolume(id1, tet_neighbours)
			var rest_vol: float = tet_volumes[id1 + tet_neighbours]
			var C: float = vol - rest_vol
			var s: float = -C / (w + alpha)
			pos[id1] += grad1 * s * inv_mass[id1]
			
			tet_neighbours += 3

func getTetVolume(base_index: int, stride: int) -> float:
	var id1: int = base_index
	var id2: int = tet_Ids[id1 + stride]
	var id3: int = tet_Ids[id1 + stride + 1]
	var id4: int = tet_Ids[id1 + stride + 2]
	
	var vec_21: Vector3 = pos[id2] - pos[id1]
	var vec_31: Vector3 = pos[id3] - pos[id1]
	var vec_41: Vector3 = pos[id4] - pos[id1]
	
	var grad: Vector3 = vec_21.cross(vec_31)
	return grad.dot(vec_41) / 6.0

func postSolve(dt):
	for i in range(pos.size()):
		if inv_mass[i] == 0.0:
			continue
		velocity[i] = (pos[i] - prev_pos[i]) * (1.0 / dt)
	# update mesh

func computeEdgeRestLengths():
	for i in range(edge_vertex_indices.size()):
		var id1: int = edge_Ids[edge_vertex_indices[i]]
		var id1_index: int = edge_vertex_indices[i]
		var s: int = edge_stride[i]
		edge_lengths[id1_index] = 0.0
		for str in range(1, s + 1):
			var id2: int = edge_Ids[id1 + str]
			var id2_index: int = id1_index + str
			var vert1: Vector3 = pos[id1]
			var vert2: Vector3 = pos[id2]
			var len: float = (vert1 - vert2).length()
			edge_lengths[id2_index] = len

func computeTetRestVolumes():
	for i in range(tet_vertex_indices.size()):
		var id1: int = tet_Ids[tet_vertex_indices[i]]
		var id1_index: int = tet_vertex_indices[i]
		var stride: float = tet_stride[i]
		tet_volumes[id1_index] = 0.0
		
		var num_iter: int = stride / 3
		var tet_neighbours: int = 1
		for iter in range(num_iter):
			var id2_index: int = id1_index + tet_neighbours
			var id3_index: int = id1_index + tet_neighbours + 1
			var id4_index: int = id1_index + tet_neighbours + 2
			
			var vol: float = getTetVolume(id1, tet_neighbours)
			
			tet_volumes[id2_index] = vol
			tet_volumes[id3_index] = vol
			tet_volumes[id4_index] = vol
			
			tet_neighbours += 3

func _ready():
	tet_mesh = TetGenMesh.new()
	tet_mesh.load_from_base_name("res://Mesh/Suzanne")
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
		edge_Ids = tet_array_constructor.construct_compact_edge(tet_mesh.edges, edge_stride, edge_vertex_indices)
		tet_Ids = tet_array_constructor.construct_compact_tet(tet_mesh.tetrahedra, tet_stride, tet_vertex_indices)
		edge_lengths.resize(edge_Ids.size())
		tet_volumes.resize(tet_Ids.size())
		tet_volumes.fill(0.0)
		computeEdgeRestLengths()
		computeTetRestVolumes()
		
		print("Edge Ids: ", edge_Ids.slice(0, 20))
		print("Tet Ids: ", tet_Ids.slice(0, 20))
		print("Edge Stride: ", edge_stride.size())
		print("Edge Vertex Indices: ", edge_vertex_indices.size())
		print("Edge Rest Lengths: ", edge_lengths.slice(0, 20))
		print("Tet Stride: ", tet_stride.size())
		print("Tet Vertex Indices: ", tet_vertex_indices.slice(0, 20))
		print("Tet Rest Volumes: ", tet_volumes.slice(0, 20))
		print("Pos: ", pos.size())
