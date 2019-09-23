using System.Threading.Tasks;

namespace BlazorServer.Core
{
    public class No_liniear_problem
    {
        MeshBuilder mesh;

        public No_liniear_problem()
        {
        }

        public void Start_calculate()
        {
        }

        public Task<MeshBuilder> GetResult()
        {
            var creator = new Gmsh3DMeshCreator();
            return Task.FromResult(creator.Initialize("Data/data.msh"));
        }
    }
}