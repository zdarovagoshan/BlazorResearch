using System.Threading.Tasks;

namespace BlazorServer.Core
{
    public class Gmsh3DMeshCreator
    {
        MeshBuilder meshBuilder;
        public int Order { get; set; } = 1;
        public MeshBuilder Initialize(string filename)
        {
            using (var stream = System.IO.File.Open(filename, System.IO.FileMode.Open, System.IO.FileAccess.Read))
            {
                meshBuilder = GmshMeshImporter.LoadMesh(stream, Order);
                return meshBuilder;
            }
        }
    }
}