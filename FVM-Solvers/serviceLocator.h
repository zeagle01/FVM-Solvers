#ifndef SERVICELOCATOR_H_H
#define SERVICELOCATOR_H_H


#include <memory>
#include <string>
#include <map>
#include <iostream>
using namespace std;


///////////////////////////////////////////////someone's implementaion of SL//////////////////////////////////////////////////
struct interface_t
{
	virtual ~interface_t() {}
};

template<class T>
class member_t
{
public:
	static interface_t* create()
	{
		return new T();
	}
};

class servicelocator_t
{
public:
	typedef map<string, interface_t*(*)()>    _classes_t;
	_classes_t                                _classes;
	typedef map<string, interface_t*>        _singletons_t;
	_singletons_t                            _singletons;
	~servicelocator_t()
	{
		for (_singletons_t::iterator it = _singletons.begin(); it != _singletons.end(); it++)
			delete it->second;
	}

	template<class T>
	void register_class(const string& id)
	{//注册
		_classes.insert(make_pair(id, T::create));
	}


	// Gets a transient instance: each call creates and returns a reference to a new object
	template<class T>
	auto_ptr<T> get_new_instance(const string& id)
	{
		_classes_t::iterator found = _classes.find(id);//在map中寻找
		if (found != _classes.end())
			return auto_ptr<T>(dynamic_cast<T*>(found->second()));
		throw runtime_error("invalid id: " + id);
	}


	template<class T>
	T* get_single_instance(const string& id)
	{//返回单例实例
		_singletons_t::iterator found_singleton = _singletons.find(id);
		//先在单例集合中寻找，
		if (found_singleton != _singletons.end())
			return dynamic_cast<T*>(found_singleton->second);
		//若单例集合中没有，则在注册集合中找,找到则放入单例集合中
		_classes_t::iterator found_class = _classes.find(id);
		if (found_class != _classes.end())
		{
			T* obj(dynamic_cast<T*>(found_class->second()));
			_singletons.insert(make_pair(id, obj));
			return obj;
		}
		throw runtime_error("invalid id: " + id);
	}

};


#endif