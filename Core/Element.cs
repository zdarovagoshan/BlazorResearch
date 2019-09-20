namespace BlazorServer.Core
{
    public class Element
    {
        public int[] Indices { get; set; }
        public int MaterialNumber { get; set; }
        public int DomainNumber { get; set; }
        public int ElementOrder { get; set; }
        public ElementType Type { get; set; }
    }

    public enum ElementType{
        Point = 0,
        Segment = 1,
        Triangle = 2,
        Quadrangle =3
    }
}