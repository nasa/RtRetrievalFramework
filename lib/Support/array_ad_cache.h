#ifndef ARRAY_AD_CACHE_H
#define ARRAY_AD_CACHE_H

#include <map>
#include "array_ad.h"

namespace FullPhysics {
  
  /****************************************************************//**
   Defines an interface for caching ArrayAd objects such that the 
   underlying caching mechanism is not necessarily exposed to
   a user of the cache. Also alternative caching mechanism
   can be switched in without a cache user needing to change anything.
  *******************************************************************/
                                                                    
  template<class K, class T, int D> class ArrayAdCache
  {
    public:
      virtual ~ArrayAdCache() {}
      virtual ArrayAd<T, D> operator[](const K& key) const = 0;
      virtual bool is_valid(const K& key) const = 0;
      virtual void insert(const K& key, const ArrayAd<T, D>& value) = 0;
      virtual void erase(const K& key) = 0;
      virtual void clear() = 0;
  };

  /****************************************************************//**
   Caches ArrayAd objects using a std::map.
   Best used when items are inserted and retrieved in an unknown order.
   *******************************************************************/
  
  template<class K, class T, int D> class ArrayAdMapCache : public ArrayAdCache<K, T, D>
  {
    public:
      ArrayAdMapCache() { }
      virtual ~ArrayAdMapCache() {}
      inline virtual ArrayAd<T, D> operator[](const K& key) const
      { 
        typename std::map<K, ArrayAd<T, D> >::const_iterator value_it;
        value_it = cache_map.find(key);
        if (value_it == cache_map.end())
          throw Exception("Could not find value_it in cache");
        return value_it->second;
      }

      inline virtual bool is_valid(const K& key) const
      { 
        return cache_map.find(key) != cache_map.end();
      }

      inline virtual void insert(const K& key, const ArrayAd<T, D>& value)
      { 
        ArrayAd<T, D> cached_copy;
        cached_copy = value;
        cache_map.insert(std::pair<K, ArrayAd<T, D> >(key, cached_copy));
      }

      inline virtual void erase(const K& key)
      {
        cache_map.erase(key);
      }

      inline virtual void clear()
      {
        cache_map.clear();
      }
    private:
      std::map<K, ArrayAd<T, D> > cache_map;
  };

  /****************************************************************//**
   Caches ArrayAd objects using a std::vector.
   This class should be used when the cached data is calculated in
   order and will be used again in order. Using operator[] is 
   slower than the same one from a map cache.
   
   The cached data is available through a method so it can be looped
   over. Which is what you want anyways.

   Additional speed can be obtained if the size of the data to cache
   is known beforehand and that size is passed during construction. The
   class will not release memory unless it is specifically instructed
   and will reuse the vector on subsequent clears and insertion cycles.
   *******************************************************************/

  template<class K, class T, int D> class ArrayAdVectorCache : public ArrayAdCache<K, T, D>
  {
    public:

      /// Optionally pass a preallocation size
      ArrayAdVectorCache(int prealloc_size = 0) : cache_index(0) 
      {
        cached_data_.resize(prealloc_size);
      }

      /// Finds the requested item using a STL lower bound call. 
      inline typename std::vector<std::pair<K, ArrayAd<T, D> > >::const_iterator find(const K& key) const
      {
        typename std::vector<std::pair<K, ArrayAd<T, D> > >::const_iterator curr_it;
        curr_it = cached_data_.begin();
        for ( ; curr_it != cached_data_.end() && curr_it != cached_data_.begin()+cache_index; curr_it++) if ( curr_it->first == key ) break;
        if(curr_it == cached_data_.begin()+cache_index)
          return cached_data_.end();
        else
          return curr_it;
      }

      /// Finds the requested item using a STL lower bound call. 
      inline typename std::vector<std::pair<K, ArrayAd<T, D> > >::iterator find(const K& key)
      {
        typename std::vector<std::pair<K, ArrayAd<T, D> > >::iterator curr_it;
        curr_it = cached_data_.begin();
        for ( ; curr_it != cached_data_.end() && curr_it != cached_data_.begin()+cache_index; curr_it++) if ( curr_it->first == key ) break;
        if(curr_it == cached_data_.begin()+cache_index)
          return cached_data_.end();
        else
          return curr_it;
      }

      /// Retrieve a specific key and error if it is no keys are available
      /// However, will return the closest one available
      inline virtual ArrayAd<T, D> operator[](const K& key) const
      {
        typename std::vector<std::pair<K, ArrayAd<T, D> > >::const_iterator value_it;
        value_it = find(key);
        if (value_it == cached_data_.end())
          throw Exception("Could not find key in vector cache to retrieve");
        return value_it->second;
      }

      inline virtual bool is_valid(const K& key) const
      { 
        return find(key) != cached_data_.end();
      }

      /// Insert the cached item into the currently available slot in the vector
      /// or creates a new slot.
      inline virtual void insert(const K& key, const ArrayAd<T, D>& value)
      { 
        ArrayAd<T, D> cached_copy;
        cached_copy = value;
        if ((int) cached_data_.size() > cache_index)
          cached_data_[cache_index] = std::pair<K, ArrayAd<T, D> >(key, cached_copy);
        else
          cached_data_.push_back(std::pair<K, ArrayAd<T, D> >(key, cached_copy));
        cache_index++;
      }

      /// Removes a specific item from the vector
      /// This call can be expensive since it will need to shift items down
      inline virtual void erase(const K& key)
      {
        typename std::vector<std::pair<K, ArrayAd<T, D> > >::iterator value_it;
        value_it = find(key);
        if (value_it == cached_data_.end())
          throw Exception("Could not find key in vector cache to erase");
        cached_data_.erase(value_it);
        cache_index--;
      }

      /// clear simply resets the location of of the cache index by default
      /// without clearing any memory.
      inline virtual void clear()
      {
        cache_index = 0;
      }

      /// Clear the vector as well as reset insertion location
      inline virtual void clear(bool clear_memory)
      {
        if (clear_memory)
          cached_data_.clear();
        clear();
      }

      /// Obtain cached data up to last point inserted.
      virtual std::vector<std::pair<K, ArrayAd<T, D> > > cached_data() const
      {
         return std::vector<std::pair<K, ArrayAd<T, D> > >(cached_data_.begin(), cached_data_.begin()+cache_index);
      }
    private:
      int cache_index;
      std::vector<std::pair<K, ArrayAd<T, D> > > cached_data_;
  };

}

#endif
