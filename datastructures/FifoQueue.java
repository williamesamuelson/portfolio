package queue_singlelinkedlist;

import java.util.*;

public class FifoQueue<E> extends AbstractQueue<E> implements Queue<E> {
	private QueueNode<E> last;
	private int size;

	public FifoQueue() {
		super();
		last = null;
		size = 0;
	}

	/**
	 * Inserts the specified element into this queue, if possible post: The
	 * specified element is added to the rear of this queue
	 * 
	 * @param e the element to insert
	 * @return true if it was possible to add the element to this queue, else false
	 */
	public boolean offer(E e) {
		QueueNode<E> qn = new QueueNode<E>(e);
		if (size == 0) {
			last = qn;
			last.next = qn;
		} else {
			QueueNode<E> first = last.next;
			last.next = qn;
			last = qn;
			qn.next = first;
		}
		size++;
		return true;
	}

	/**
	 * Returns the number of elements in this queue
	 * 
	 * @return the number of elements in this queue
	 */
	public int size() {
		return size;
	}

	/**
	 * Retrieves, but does not remove, the head of this queue, returning null if
	 * this queue is empty
	 * 
	 * @return the head element of this queue, or null if this queue is empty
	 */
	public E peek() {
		if (size == 0) {
			return null;
		} else {
			return last.next.element;
		}
	}

	/**
	 * Retrieves and removes the head of this queue, or null if this queue is empty.
	 * post: the head of the queue is removed if it was not empty
	 * 
	 * @return the head of this queue, or null if the queue is empty
	 */
	public E poll() {
		if (size == 0) {
			return null;
		} else if (size == 1) {
			E temp = last.element;
			last = null;
			size--;
			return temp;
		} else {
			QueueNode<E> first = last.next;
			last.next = first.next;
			size--;
			return first.element;
		}
	}
	
	/**
	* Appends the specified queue to this queue
	* post: all elements from the specified queue are appended
	* to this queue. The specified queue (q) is empty after the call.
	* @param q the queue to append
	* @throws IllegalArgumentException if this queue and q are identical
	*/
	public void append(FifoQueue<E> q) {
		if (q == this) {
			throw new IllegalArgumentException();
		} else if (q.size != 0 && this.size == 0) {
			last = q.last;
			this.size = q.size;
			q.size = 0;
			q.last = null;
		} else if (this.size != 0 && q.size != 0) {
			QueueNode<E> first = last.next;
			last.next = q.last.next;
			last = q.last;
			last.next = first;
			this.size += q.size;
			q.size = 0;
			q.last = null;
		}
	}

	
	
	/**
	 * Returns an iterator over the elements in this queue
	 * 
	 * @return an iterator over the elements in this queue
	 */
	public Iterator<E> iterator() {
		return new QueueIterator();
	}

	private class QueueIterator implements Iterator<E> {
		private QueueNode<E> pos;
		private boolean atStart;

		private QueueIterator() {
			if (size == 0) {
				pos = null;
			} else {
				pos = last.next;
				atStart = true;
			}
		}

		public boolean hasNext() {
			if (size == 0 || pos == last.next && !atStart) {
				return false;
			}
			return true;
		}

		public E next() {
			if (!hasNext()) {
				throw new NoSuchElementException();
			} else {
				E temp = pos.element;
				pos = pos.next;
				atStart = false;
				return temp;
			}
		}
	}

	private static class QueueNode<E> {
		E element;
		QueueNode<E> next;

		private QueueNode(E x) {
			element = x;
			next = null;
		}
	}

}
