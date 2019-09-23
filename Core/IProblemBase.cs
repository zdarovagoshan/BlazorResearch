using System.Collections.Generic;

namespace BlazorServer.Core
{
    public interface IClearable
    {
        void Clear();
    }
    public interface IMaterialCatalog<T> : IClearable where T : IMaterialDescriptor
    {
        bool TryGetValue(int cat, out T value);
        bool TryGetValue(string name, out T value);
        void Add(T desc);
        void AddDefault(string name);
        void Rename(string oldValue, string newValue);
        void Renumber(int oldValue, int newValue);
        bool Remove(string name);
        IEnumerable<T> Enumerate { get; }
    }
    public interface ICatalogManager : IClearable
    {
        // IMaterialCatalog<IVolumeMaterialDescriptor> MatCat { get; }
        // IMaterialCatalog<IBoundaryMaterialDescriptor> BoundCat { get; }
        int AddMaterial(string name);
        int AddBoundary(string name);
    }
}