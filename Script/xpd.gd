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
		var id1: int = edge_vertex_indices[i]
		var stride = edge_stride[i]
		for stri in range(1, stride + 1):
			var id2: int = edge_vertex_indices[id1 + stri]
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
# [3, 6] # Tetrahedron Stride
# [0, 4] # Tetrahedron Vertex Indices
func solveVolumes(compliance: float, dt: float):
	pass

func postSolve(dt):
	for i in range(pos.size()):
		if inv_mass[i] == 0.0:
			continue
		velocity[i] = (pos[i] - prev_pos[i]) * (1.0 / dt)
	# update mesh

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
		
		print("Edge Stride: ", edge_stride.size())
		print("Edge Vertex Indices: ", edge_vertex_indices.size())
		print("Tet Stride: ", tet_stride.size())
		print("Tet Vertex Indices: ", tet_vertex_indices.size())
