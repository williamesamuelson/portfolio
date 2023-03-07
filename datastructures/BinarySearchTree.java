package bst;

public class BinarySearchTree<E extends Comparable<? super E>> {
	BinaryNode<E> root;
	int size;

	public static void main(String[] args) {
		BinarySearchTree<Integer> bst = new BinarySearchTree<Integer>();
		for (int i = 0; i < 5; i++) {
			bst.add(i);
		}
		for (int i = -5; i < 0; i++) {
			bst.add(i);
		}
		bst.rebuild();
		BSTVisualizer vis = new BSTVisualizer("Binary Search Tree", 400, 200);
		vis.drawTree(bst);
	}

	/**
	 * Constructs an empty binary searchtree.
	 */
	public BinarySearchTree() {
	}

	/**
	 * Inserts the specified element in the tree if no duplicate exists.
	 * 
	 * @param x element to be inserted
	 * @return true if the the element was inserted
	 */
	public boolean add(E x) {
		return add(x, root);
	}

	private boolean add(E x, BinaryNode<E> n) {
		if (root == null) {
			root = new BinaryNode<E>(x);
			size++;
			return true;
		} else if (n.right == null && x.compareTo(n.element) > 0) {
			n.right = new BinaryNode<>(x);
			size++;
			return true;
		} else if (n.left == null && x.compareTo(n.element) < 0) {
			n.left = new BinaryNode<>(x);
			size++;
			return true;
		} else {
			if (x.compareTo(n.element) > 0) {
				return add(x, n.right);
			} else if (x.compareTo(n.element) < 0) {
				return add(x, n.left);
			} else {
				return false;
			}
		}
	}

	/**
	 * Computes the height of tree.
	 * 
	 * @return the height of the tree
	 */
	public int height() {
		return height(root);
	}

	private int height(BinaryNode<E> n) {
		if (n == null) {
			return 0;
		} else {
			return 1 + Math.max(height(n.left), height(n.right));
		}
	}

	/**
	 * Returns the number of elements in this tree.
	 * 
	 * @return the number of elements in this tree
	 */
	public int size() {
		return size;
	}

	/**
	 * Print tree contents in inorder.
	 */
	public void printTree() {
		printTree(root);
	}

	private void printTree(BinaryNode<E> n) {
		if (n != null) {
			printTree(n.left);
			System.out.println(n.element);
			printTree(n.right);
		}
	}

	/**
	 * Builds a complete tree from the elements in the tree.
	 */
	public void rebuild() {
		@SuppressWarnings("unchecked")
		E[] a = (E[]) new Comparable[size];
		toArray(root, a, 0);
		root = buildTree(a, 0, a.length - 1);
	}

	/*
	 * Adds all elements from the tree rooted at n in inorder to the array a
	 * starting at a[index]. Returns the index of the last inserted element + 1 (the
	 * first empty position in a).
	 */
	private int toArray(BinaryNode<E> n, E[] a, int index) {
		if (index < a.length) {
			if (n.left != null && n.right != null) {
				int q = toArray(n.left, a, index);
				a[q] = n.element;
				int p = toArray(n.right, a, q + 1);
				return p;
			} else if (n.left != null && n.right == null) {
				int q = toArray(n.left, a, index);
				a[q] = n.element;
				return q + 1;
			} else if (n.left == null && n.right != null) {
				a[index] = n.element;
				int q = toArray(n.right, a, index + 1);
				return q;
			} else {
				a[index] = n.element;
				return index + 1;
			}
		} else {
			return a.length;
		}
	}

	/*
	 * Builds a complete tree from the elements a[first]..a[last]. Elements in the
	 * array a are assumed to be in ascending order. Returns the root of tree.
	 */
	private BinaryNode<E> buildTree(E[] a, int first, int last) {
		if (first > last) {
			return null;
		} else {
			int mid = first + ((last - first) / 2);
			BinaryNode<E> bn = new BinaryNode<E>(a[mid]);
			BinaryNode<E> bnLeft = buildTree(a, first, mid - 1);
			BinaryNode<E> bnRight = buildTree(a, mid + 1, last);
			bn.left = bnLeft;
			bn.right = bnRight;
			return bn;
		}
	}

	static class BinaryNode<E> {
		E element;
		BinaryNode<E> left;
		BinaryNode<E> right;

		private BinaryNode(E element) {
			this.element = element;
		}
	}

}
