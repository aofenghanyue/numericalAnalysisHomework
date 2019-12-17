#pragma once
#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdio.h> // NULL

using namespace std;

template <typename T>
class Vector {/*�����࣬Ϊ�ӿ��ٶȣ����ɶ�̬�����ڴ�*/
private:
	int length;
	T* elements;

public:
	Vector();
	explicit Vector(int n);
	Vector(int n, const T& a);
	Vector(int n, const T* a);
	Vector(const Vector& other_vector);
	Vector& operator=(const Vector& other_vector);
	inline T& operator[](const int i);
	inline const T& operator[](const int i) const;
	inline int size() const;
	inline Vector<T> operator/(T b);
	inline Vector<T> operator*(T b);
	~Vector();
};

template<typename T>
Vector<T>::Vector() : length(0), elements(NULL) {
	/*Ĭ�Ϲ�������Ϊ0��Ԫ��Ϊ�յ�����*/
}

template<typename T>
Vector<T>::Vector(int n) : length(n), elements(n > 0 ? new T[n] : NULL) {
	/*��������Ϊn��Ԫ��Ϊ�յ�����*/
}

template<typename T>
Vector<T>::Vector(int n, const T& a) : length(n), elements(n > 0 ? new T[n] : NULL) {
	/*��������Ϊn��Ԫ��Ϊa������*/
	for (int i = 0; i < n; i++) elements[i] = a;
}

template<typename T>
Vector<T>::Vector(int n, const T* a) : length(n), elements(n > 0 ? new T[n] : NULL) {
	/*��������Ϊn��Ԫ�غ�a��ָ��������ͬ������*/
	for (int i = 0; i < n; i++) elements[i] = *a++;
}

template<typename T>
Vector<T>::Vector(const Vector& other_vector) : length(other_vector.length), elements(length > 0 ? new T[length] : NULL) {
	/*����һ������*/
	for (int i = 0; i < length; i++) elements[i] = other_vector[i];
}

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other_vector) {
	/*����һ������
	* ͨ�������� = ��ʵ��
	* ����һ���µĶ���
	* ���磺Vector v = v0;
	*/
	if (this->length != other_vector.length) {
		if (elements != NULL) delete[](elements);
		length = other_vector.length;
		elements = new T[length];
	}
	for (int i = 0; i < length; i++) {
		elements[i] = other_vector[i];
	}

	return *this;
}

template<typename T>
inline T& Vector<T>::operator[](const int i) {
	return elements[i];
}

template<typename T>
inline const T& Vector<T>::operator[](const int i) const {
	return elements[i];
}

template<typename T>
inline int Vector<T>::size() const {
	/*����������С*/
	return length > 0 ? length : 0;
}


template <typename T>
inline Vector<T> Vector<T>::operator/(T b) {
	Vector <T> result(length);
	for (int i = 0; i < length; i++) {
		result[i] = elements[i] / b;
	}
	return result;
}

template <typename T>
inline Vector<T> Vector<T>::operator*(T b) {
	Vector <T> result(length);
	for (int i = 0; i < length; i++) {
		result[i] = elements[i] * b;
	}
	return result;
}

template<typename T>
inline Vector<T>::~Vector() {
	if (elements) delete[] elements;
}

typedef Vector<double> VectorD;

#endif