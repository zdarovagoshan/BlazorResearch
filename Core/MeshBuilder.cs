using System;
using System.Collections.Generic;
using System.Linq;

namespace BlazorServer.Core
{
    public class MeshBuilder
    {
        static readonly Element Zero = new Element { Type = (ElementType)(-1), MaterialNumber = -1, DomainNumber = -1 };
        SimpleEnumerator<Vector3D> _vertices = new SimpleEnumerator<Vector3D>();
        List<Element> _elements = new List<Element>();
        public int _externMaterial = -1;
        public Func<IEnumerable<Vector3D>, int, int, int> BoundaryMaterialSelector { get; set; }
        public IReadOnlyList<Vector3D> Vertices { get; set; }
        public IEnumerable<Element> Elements => _elements.Skip(1);
        public int Dimension { get; set; }

        public IReadOnlyDictionary<int, string> VolumeMaterialNames { get; set; }
        public IReadOnlyDictionary<int, string> BoundaryMaterialNames { get; set; }
        public MeshBuilder()
        {
            _elements.Add(Zero);
        }

        public void AddElements(IEnumerable<Element> elements, bool repair = true)
        {
            foreach (var elem in elements)
            {
                AddElement(elem, repair);
            }
        }
        public void AddElement(Element element, bool repair = true)
        {
            if (element.Type != ElementType.Point)
            {
                if (repair)
                {
                    RepairElement(element);
                    _elements.Add(element);
                }
            }
        }
        private void RepairElement(Element element)
        {
            var indices = element.Indices;
            var pts = indices.Select(i => Vertices[i]).ToArray();
            switch (element.Type)
            {
                case ElementType.Quadrangle:
                    {
                        Vector3D n1 = Vector3D.Cross(pts[1] - pts[0], pts[3] - pts[1]);
                        Vector3D n2 = Vector3D.Cross(pts[3] - pts[0], pts[2] - pts[3]);
                        if (n1 * n2 < 0)
                            (indices[2], indices[3]) = (indices[3], indices[2]);
                    }
                    break;

                default:
                    break;
            }
        }

    }
}