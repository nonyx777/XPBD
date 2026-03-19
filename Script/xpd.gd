extends Node

var data = JSON.parse_string(FileAccess.get_file_as_string("res://tetrahedral_mesh.json"))

var vertices = data["vertices"]
var tetIds = data["tetIds"]
var tetEdgeIds = data["tetEdgeIds"]
var surface = data["tetSurfaceTriIds"]

func _ready() -> void:
	var verts = PackedVector3Array()
	for v in vertices:
		verts.append(Vector3(v[0], v[1], v[2]))
	var indices = PackedInt32Array()
	for tri in surface:
		indices.append(tri[0])
		indices.append(tri[1])
		indices.append(tri[2])
	
	var normals = PackedVector3Array()
	normals.resize(verts.size())
	for i in range(normals.size()):
		normals[i] = Vector3.ZERO
	for i in range(0, indices.size(), 3):
		var i0 = indices[i]
		var i1 = indices[i+1]
		var i2 = indices[i+2]
		
		var n = (verts[i1] - verts[i0]).cross(verts[i2] - verts[i0])
		normals[i0] += n
		normals[i1] += n
		normals[i2] += n
	for i in range(normals.size()):
		normals[i] = normals[i].normalized()
	
	var arrays = []
	arrays.resize(Mesh.ARRAY_MAX)
	arrays[Mesh.ARRAY_VERTEX] = verts
	arrays[Mesh.ARRAY_INDEX] = indices
	arrays[Mesh.ARRAY_NORMAL] = normals
	
	var mesh = ArrayMesh.new()
	mesh.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	
	var mi = MeshInstance3D.new()
	mi.mesh = mesh
	add_child(mi)
	
	print("Vert Size: ", verts.size())
	print("Indices Size: ", indices.size())
	print("Normals Size: ", normals.size())
