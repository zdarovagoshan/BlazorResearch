using System.Collections;
using System.Collections.Generic;

namespace BlazorServer.Core
{
    public class SimpleEnumerator<T> : IReadOnlyList<T>
    {
        List<T> list = new List<T>();
        Dictionary<T, int> dic;
        public SimpleEnumerator()
            : this(EqualityComparer<T>.Default)
        {
        }
        public SimpleEnumerator(IEqualityComparer<T> comparer)
        {
            dic = new Dictionary<T, int>(comparer);
        }
        public int IndexOf(T item)
        {
            int res;
            if (!dic.TryGetValue(item, out res))
            {
                dic[item] = res = list.Count;
                list.Add(item);
            }
            return res;
        }
        public T this[int index] { get { return list[index]; } }
        public int Count { get { return list.Count; } }
        public IEnumerator<T> GetEnumerator()
        {
            return list.GetEnumerator();
        }
        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}