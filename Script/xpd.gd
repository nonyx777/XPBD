extends Node3D

var tet_mesh: TetGenMesh
var tet_array_constructor: TetArrayConstructor

var verts: PackedVector3Array
var edgeIds: PackedInt32Array
var tetIds: PackedInt32Array

var edge_stride: PackedInt32Array
var edge_vertex_indices: PackedInt32Array

var tet_stride: PackedInt32Array
var tet_vertex_indices: PackedInt32Array

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
		
		verts = tet_mesh.vertices
		edgeIds = tet_array_constructor.construct_compact_edge(tet_mesh.edges, edge_stride, edge_vertex_indices)
		tetIds = tet_array_constructor.construct_compact_tet(tet_mesh.tetrahedra, tet_stride, tet_vertex_indices)
		
		print("Edge Stride: ", edge_stride.size())
		print("Edge Vertex Indices: ", edge_vertex_indices.size())
		print("Tet Stride: ", tet_stride.size())
		print("Tet Vertex Indices: ", tet_vertex_indices.size())
