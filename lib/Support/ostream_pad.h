#ifndef OSTREAM_PAD_H
#define OSTREAM_PAD_H
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/foreach.hpp>

namespace FullPhysics {
class OstreamPadFilter : public boost::iostreams::output_filter {
public:
  OstreamPadFilter(const std::string& Pad) 
  : pad(Pad), need_pad(true) {}
  template <typename Sink>
  bool put(Sink& snk, char c)
  {
    if(need_pad) {
      BOOST_FOREACH(char c2, pad)
	boost::iostreams::put(snk, c2);
    }
    need_pad = (c == '\n');
    return boost::iostreams::put(snk, c);
  }
private:
  std::string pad;
  bool need_pad;
};

/****************************************************************//**
  This is a filtering stream that adds a pad to the front of every
  line written out. This can be used to do simple formating (among
  other things) of adding space for a nested item.
*******************************************************************/

class OstreamPad : public boost::iostreams::filtering_ostream {
public:
  OstreamPad(std::ostream& Os, const std::string& Pad)
    : os(Os), p(Pad) { push(p); push(Os); }
private:
  std::ostream& os;
  OstreamPadFilter p;
};
}
#endif
