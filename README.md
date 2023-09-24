# Binary Heap Implementation

This repository hosts a binary heap implementation, an essential data structure used for priority queues. A binary heap is a complete binary tree with the remarkable feature that the value stored at each node is less than or equal to the values of its children, ensuring efficient retrieval of elements with the highest or lowest priority.

## Features

This binary heap implementation offers the following features:

- *Shape Property:* The binary heap is a complete binary tree. All levels are filled, except possibly the last level, which is filled from left to right.

- *Heap Property:* The value at each node is less than or equal to the values of its children, maintaining a heap property that allows for efficient prioritization.

- *Operations:* The following operations are supported:
  - `CPBH IsEmpty`: Check if the binary heap is empty.
  - `CPBH Enqueue`: Add a new node while preserving the shape property.
  - `CPBH ReEnqueue`: Update the value of an existing node while maintaining the heap property.
  - `CPBH Dequeue`: Retrieve and remove the root node, maintaining the heap's integrity.

## Data Representation

This implementation uses a memory-efficient data representation that leverages an array of arrays to manage each level of the binary heap. This structure allows for practical indexing and efficient management of heaps of varying sizes.

## Usage

This binary heap implementation can be used in a variety of applications that require priority queue functionality. It provides a solid foundation for tasks that involve ordering elements by priority, such as in graph algorithms, scheduling, and more.
