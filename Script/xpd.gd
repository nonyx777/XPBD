extends Node3D

var tet_mesh: TetGenMesh

func _ready():
	tet_mesh = TetGenMesh.new()
	tet_mesh.load_from_base_name("res://Mesh/Suzanne")
	
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
