package map;

public class SimpleHashMap<K, V> implements Map<K, V> {
	private Entry<K, V>[] table;
	private int size;
	private int capacity;
	private double loadFactor;

	@SuppressWarnings("unchecked")
	public SimpleHashMap() {
		this(16);
	}

	@SuppressWarnings("unchecked")
	public SimpleHashMap(int capacity) {
		this.capacity = capacity;
		loadFactor = 0.75;
		size = 0;
		table = (Entry<K, V>[]) new Entry[capacity];
	}

	@Override
	public V get(Object object) {
		@SuppressWarnings("unchecked")
		K key = (K) object;
		Entry<K, V> e = find(index(key), key);
		if (e != null) {
			return e.value;
		} else {
			return null;
		}
	}

	@Override
	public boolean isEmpty() {
		return size == 0;
	}

	@Override
	public V put(K key, V value) {
		if ((double) (size + 1) / capacity > loadFactor) {
			rehash();
		}
		int index = index(key);
		Entry<K, V> entry = find(index, key);
		if (entry != null) {
			V temp = entry.value;
			entry.setValue(value);
			return temp;
		} else {
			Entry<K, V> first = table[index];
			Entry<K, V> newEntry = new Entry<K, V>(key, value);
			table[index] = newEntry;
			size++;
			newEntry.next = first;
			return null;
		}
	}

	@Override
	public V remove(Object obj) {
		@SuppressWarnings("unchecked")
		K key = (K) obj;
		int index = index(key);
		Entry<K, V> first = table[index];
		if (first == null) {
			return null;
		} else if (first.key.equals(key)) {
			table[index] = first.next;
			size--;
			return first.value;
		} else {
			Entry<K, V> current = first.next;
			Entry<K, V> pre = first;
			while (current != null) {
				if (current.key.equals(key)) {
					pre.next = current.next;
					size--;
					return current.value;
				}
				pre = current;
				current = current.next;
			}
			return null;
		}

	}

	@Override
	public int size() {
		return size;
	}

	public String show() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < capacity; i++) {
			sb.append(i + "    ");
			Entry<K, V> current = table[i];
			while (current != null) {
				sb.append(current.toString() + " ");
				current = current.next;
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	private int index(K key) {
		return Math.abs(key.hashCode() % capacity);
	}

	private Entry<K, V> find(int index, K key) {
		Entry<K, V> q = table[index];
		while (q != null) {
			if (q.key.equals(key)) {
				return q;
			}
			q = q.next;
		}
		return null;
	}

	private void rehash() {
		@SuppressWarnings("unchecked")
		Entry<K, V>[] newTable = (Entry<K, V>[]) new Entry[capacity * 2];
		int oldCap = capacity;
		capacity *= 2;
		for (int i = 0; i < oldCap; i++) {
			Entry<K, V> current = table[i];
			while (current != null) {
				int newIndex = index(current.key);
				Entry<K, V> next = current.next;
				Entry<K, V> first = newTable[newIndex]; //first == nuvarande första i nya tabellen
				current.next = first;                   
				newTable[newIndex] = current;           //byter ut first mot nya entryn
				current = next;

			}
		}
		table = newTable;
	}

	public static class Entry<K, V> implements Map.Entry<K, V> {
		private K key;
		private V value;
		private Entry<K, V> next;

		public Entry(K key, V value) {
			this.key = key;
			this.value = value;
		}

		@Override
		public K getKey() {
			return key;
		}

		@Override
		public V getValue() {
			return value;
		}

		@Override
		public V setValue(V value) {
			V temp = value;
			this.value = value;
			return temp;
		}

		public String toString() {
			return key + "=" + value;
		}
	}

}
