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

func solveEdges(compliance: float, dt: float):
	pass
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
