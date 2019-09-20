using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace BlazorServer.Core
{
    public class GmshMeshImporter
    {
        enum GmshElementType
        {
            Point = 1,
            Line = 2,
            Triangle = 3,
            Quadrangle = 4,
            Line3 = 5,
            Triangle6 = 6

        }

        private static Vector3D[] ReadNodes4(StreamReader reader)
        {
            var line = reader.ReadLine().Split(' ');
            int numEntityBlocks = int.Parse(line[0]);
            int numNodes = int.Parse(line[1]);
            var vertices = new Vector3D[numNodes];
            for (int i = 0; i < numEntityBlocks; i++)
            {
                line = reader.ReadLine().Split(' ');
                int numNodesInBlock = int.Parse(line[3]);
                var numbers = new long[numNodesInBlock];
                for (int j = 0; j < numNodesInBlock; j++)
                    numbers[j] = long.Parse(reader.ReadLine());
                for (int j = 0; j < numNodesInBlock; j++)
                    vertices[numbers[j] - 1] = Vector3D.Parse(reader.ReadLine());
            }
            return vertices;
        }
        private static Dictionary<(int dim, int tag), int> ReadEntities(StreamReader reader)
        {
            var line = reader.ReadLine().Split(' ');
            int numPoints = int.Parse(line[0]);
            int numCurves = int.Parse(line[1]);
            int numSurfaces = int.Parse(line[2]);
            int numVolumes = int.Parse(line[3]);
            var dict = new Dictionary<(int dim, int tag), int>();
            void ParseEntity(int dim)
            {
                line = reader.ReadLine().Split(' ');
                int tag = int.Parse(line[0]);
                int addr = 7;
                if (dim == 0) addr = 4;
                int numPhysicals = int.Parse(line[addr]);
                if (numPhysicals > 1)
                    throw new Exception($"Gmsh 4 format error: numPhysicals = {numPhysicals} > 1");
                else if (numPhysicals == 1)
                    dict.Add((dim, tag), int.Parse(line[addr + 1]));
            }
            for (int i = 0; i < numPoints; ++i)
                ParseEntity(0);
            for (int i = 0; i < numCurves; ++i)
                ParseEntity(1);
            for (int i = 0; i < numSurfaces; ++i)
                ParseEntity(2);
            for (int i = 0; i < numVolumes; ++i)
                ParseEntity(3);
            return dict;
        }
        private static IEnumerable<Element> ReadElements4(StreamReader reader, Dictionary<(int dim, int tag), int> entities, int order)
        {
            var line = reader.ReadLine().Split(' ');
            int numEntityBlocks = int.Parse(line[0]);
            for (int i = 0; i < numEntityBlocks; ++i)
            {
                line = reader.ReadLine().Split(' ');
                int adder = 1;
                var tagEntity = int.Parse(line[(adder + 0) % 2]);
                var dimEntity = int.Parse(line[(adder + 1) % 2]);
                var typeEle = (GmshElementType)int.Parse(line[2]);
                var numElementsInBlock = int.Parse(line[3]);
                for (int j = 0; j < numElementsInBlock; ++j)
                {
                    var nodes = reader.ReadLine().TrimEnd(' ').Split(' ').Skip(1).Select(v => int.Parse(v) - 1).ToArray();
                    var elem = ConvertElement(typeEle, nodes);
                    if (elem != null)
                    {
                        elem.ElementOrder = order;
                        elem.DomainNumber = tagEntity;
                        elem.MaterialNumber = entities[(dimEntity, tagEntity)];
                        yield return elem;
                    }
                }

            }
        }
        private static Element ConvertElement(GmshElementType type, int[] nodes)
        {
            switch (type)
            {
                case GmshElementType.Line:
                case GmshElementType.Line3:
                    return new Element { Type = ElementType.Segment, Indices = nodes };
                case GmshElementType.Triangle:
                case GmshElementType.Triangle6:
                    return new Element { Type = ElementType.Triangle, Indices = nodes };
                case GmshElementType.Quadrangle:
                    return new Element { Type = ElementType.Quadrangle, Indices = new[] { nodes[0], nodes[1], nodes[3], nodes[2] } };
                default:
                    return null;
            }
        }
        private static IReadOnlyDictionary<int, string> ReadPhysicalNames(StreamReader reader)
        {
            int namesNumber = int.Parse(reader.ReadLine());
            var names = new Dictionary<int, string>();
            for (int i = 0; i < namesNumber; i++)
            {
                var line = reader.ReadLine().Split(new[] { ' ' }, 3);
                int index = int.Parse(line[1]);
                names[int.Parse(line[1])] = line[2].Trim('"');
            }
            return names;
        }
        public static MeshBuilder LoadMesh(Stream stream, int order)
        {
            using (var reader = new StreamReader(stream))
            {
                var builder = new MeshBuilder();
                IReadOnlyDictionary<int, string> MaterialNames = null;
                Dictionary<(int dim, int tag), int> entities = null;
                double version = 0;
                while (!reader.EndOfStream)
                {
                    var section = reader.ReadLine();
                    switch (section)
                    {
                        case "$MeshFormat":
                            var line = reader.ReadLine().Split(' ');
                            version = double.Parse(line[0]);
                            break;
                        case "$Entities":
                            entities = ReadEntities(reader);
                            break;
                        case "$Nodes":
                            builder.Vertices = ReadNodes4(reader);
                            break;
                        case "$Elements":
                            builder.AddElements(ReadElements4(reader, entities, order), true);
                            break;
                        case "$PhysicalNames":
                            MaterialNames = ReadPhysicalNames(reader);
                            break;
                        default:
                            break;
                    }

                    while (reader.Peek() != '$')
                        reader.ReadLine();
                    if (reader.ReadLine() != string.Concat("$End", section.Substring(1)))
                        throw new InvalidDataException("Mesh file corrupted!");
                }
                var Is3DMesh = builder.Elements.Any(e => e.Type.Dimension() == 3);
                Dictionary<int, string> VolumeDict = new Dictionary<int, string>();
                Dictionary<int, string> BoundaryDict = new Dictionary<int, string>();
                foreach (var elem in builder.Elements)
                {
                    if (MaterialNames.TryGetValue(elem.MaterialNumber, out var name))
                    {
                        if (elem.Type.Dimension() == 2)
                            VolumeDict[elem.MaterialNumber] = name;
                        else
                            BoundaryDict[elem.MaterialNumber] = name;
                    }
                }
                builder.BoundaryMaterialNames = BoundaryDict;
                builder.VolumeMaterialNames = VolumeDict;
                return builder;

            }
        }
    }
}