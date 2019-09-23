using System;

namespace BlazorServer.Core
{
    public interface IMaterialProperty
    {
        string Name { get; }
    }
    public interface IMaterialPropertyCollection
    {
        bool TryGetProperty(string name, out IMaterialProperty property);
    }
    public interface IMaterialDescriptor : ICloneable
    {
        int CatNum { get; set; }
        string Name { get; }
        bool IsReserved { get; }
        IMaterialPropertyCollection Coefs { get; }
    }
}