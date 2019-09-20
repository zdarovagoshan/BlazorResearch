namespace BlazorServer.Core
{
    public static class ElementTypesExtensions
    {
        private static int[] elemDimensions = { 0, 1, 2, 2, 3, 3, 3, 3 };
        public static int Dimension(this ElementType elem)
        {
            return elemDimensions[(int)elem];
        }
    }
}