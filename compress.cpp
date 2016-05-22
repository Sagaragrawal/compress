#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <queue>

typedef ptrdiff_t Int;
typedef size_t Uint;
typedef unsigned char Byte;

using std::cout;
using std::endl;

template<typename T>
struct Img
{
  Int sz0;
  Int sz1;
  std::shared_ptr<T> ptr;

  Img(Int sz0, Int sz1)
    : sz0(sz0), sz1(sz1), ptr(new T[sz0 * sz1], std::default_delete<T[]>())
  {
    std::fill_n(ptr.get(), sz0 * sz1, T());
  }

  T& operator()(Int i0, Int i1) { return ptr.get()[i0 * sz1 + i1]; }
  const T& operator()(Int i0, Int i1) const { return ptr.get()[i0 * sz1 + i1]; }
};

static const Int bytesize = 8;

struct BitView
{
  BitView() = default;
  Byte* ptr = nullptr;
  Int offset = 0;
  Int len = 0;
};

Int read_bit(BitView& bits)
{
  //cout << "(rb " << Int(bits.ptr[0]) << " " << bits.offset << ")";
  bool r = (bits.ptr[0] >> bits.offset) & 1;
  bits.offset++;
  bits.len--;
  if(bits.offset == bytesize)
  {
    bits.offset = 0;
    bits.ptr++;
  }

  return r;
}

Int read_bits(BitView& bits, Int n)
{
  Int r = 0;
  for(Int i = 0; i < n; i++) r |= (read_bit(bits) << i);
  return r;
}

struct BitArray
{
  BitArray() = default;
  std::vector<Byte> v = {0};
  Int offset = 0;

  BitView get_view()
  {
    return { &v[0], 0, Int((v.size() - 1) * bytesize + offset) };
  }

  void push(Int val)
  {
    v.back() |= (val << offset);
    offset++;
    if(offset == bytesize)
    {
      offset = 0;
      v.push_back(0);
    }
  }

  void push(BitView view) { while(view.len > 0) push(read_bit(view)); }

  void pop()
  {
    if(offset == 0)
    {
      offset = bytesize;
      v.pop_back();
    }

    offset--;
    v.back() &= ~(1 << offset);
  }
};

//TODO: comments
Img<Byte> load_pgm(const char* name)
{
  std::ifstream f(name, std::ios_base::binary);
  std::string line;
  std::getline(f, line);
  Int sz0, sz1, maxval;
  f >> sz1;
  f >> sz0;
  f >> maxval;
  auto r = Img<Byte>(sz0, sz1);
  for(Int i = 0; i < sz0 * sz1; i++)
  {
    Int val;
    f >> val;
    r.ptr.get()[i] = val;
  }

  return r;
}

void save_pgm(const Img<Byte>& img, const char* name)
{
  std::ofstream f(name, std::ios_base::binary);
  f << "P2" << endl;
  f << img.sz1 << " " << img.sz0 << endl;
  f << 255 << endl;

  for(Int i0 = 0; i0 < img.sz0; i0++)
  {
    const Byte* p = &img(i0, 0);
    for(Int i1 = 0; i1 < img.sz1 - 1; i1++) f << Int(p[i1]) << " ";
    f << Int(p[img.sz1 - 1]) << endl;
  }
}

#define for_img(img, i0, i1) \
  for(Int i0 = 0; i0 < (img).sz0; i0++) \
    for(Int i1 = 0; i1 < (img).sz1; i1++)

void runlength_encode(const Byte* src, Int n, Int max_len, std::vector<Byte>& dst)
{
  Byte val = src[0];
  Int idx = 0;
  for(Int i = 1; i < n; i++)
    if(src[i] != val || i - idx == max_len)
    {
      dst.push_back(val);
      dst.push_back(i - idx);
      val = src[i]; 
      idx = i;
    }

  dst.push_back(val);
  dst.push_back(n - idx);
}

void runlength_decode(const Byte* src, ptrdiff_t n, std::vector<Byte>& dst)
{
  for(Int i = 0; i < n; i += 2)
    std::fill_n(std::back_inserter(dst), src[i + 1], src[i]);
}

namespace huffman
{
  void create_encode_tree(Int* freq, Int freq_len, std::vector<Int>& dst)
  {
    std::priority_queue<
      std::pair<Int, Int>,
      std::vector<std::pair<Int, Int>>,
      std::greater<>> q;

    dst.clear();
    for(Int i = 0; i < freq_len; i++)
    {
      dst.push_back(-1);
      q.emplace(freq[i], i);
    }

    while(q.size() > 1)
    {
      auto a = q.top(); q.pop();
      auto b = q.top(); q.pop();
      Int new_idx = dst.size();
      dst.push_back(-1);
      dst[a.second] = new_idx;
      dst[b.second] = new_idx;
      q.emplace(a.first + b.first, new_idx);
    }
  }

  Int encoded_len(Int* freq, Int freq_len, const std::vector<Int>& tree)
  {
    Int r = 0;
    for(Int i = 0; i < freq_len; i++)
      for(Int j = i; j != -1; j = tree[j]) r += freq[i];

    return r;
  }

  struct Node
  {
    Int value;
    std::unique_ptr<Node> children[2];
    Node(Int value) : value(value) { }
    Node(Node* child0, Node* child1)
    {
      children[0].reset(child0);
      children[1].reset(child1);
    }
  };

  void print_node(Node* n)
  {
    if(n->children[0])
    {
      cout << "(";
      print_node(n->children[0].get());
      cout << " ";
      print_node(n->children[1].get());
      cout << ")";
    }
    else
      cout << n->value;
  }

  std::unique_ptr<Node> create_child_tree(const std::vector<Int>& parent_tree)
  {
    Int n = (parent_tree.size() + 1) / 2;
    std::vector<Node*> nodes(parent_tree.size() - n, nullptr);
    std::unique_ptr<Node> r;
    for(Int i = 0; i < n; i++)
    {
      auto node = new Node(i);
      for(Int j = parent_tree[i]; ; j = parent_tree[j])
        if(j == -1)
        {
          r.reset(node);
          break;
        }
        else if(nodes[j - n])
        {
          nodes[j - n]->children[1].reset(node);
          break;
        }
        else
        {
          node = new Node(node, nullptr);
          nodes[j - n] = node;
        }
    }

    return std::move(r);
  }

  struct Encoder
  {
    std::vector<BitArray> table;
    void recurse(BitArray& bits, Node* node)
    {
      if(node->children[0])
        for(Int i = 0; i < 2; i++)
        {
          bits.push(i);
          recurse(bits, node->children[i].get());
          bits.pop();
        }
      else table[node->value] = bits;
    }

    Encoder(const std::vector<Int>& encode_tree)
    {
      Int n = (encode_tree.size() + 1) / 2;
      BitArray empty;
      table.resize(n);
      recurse(empty, create_child_tree(encode_tree).get());
    }

    BitView operator()(Int val) { return table[val].get_view(); }
  };

  struct Decoder
  {
    std::unique_ptr<Node> tree;
    Decoder(const std::vector<Int>& parent_tree)
      : tree(create_child_tree(parent_tree)) { }

    Int operator()(BitView& bits)
    {
      Node* n = tree.get();
      while(true)
        if(n->children[0])
          n = n->children[read_bit(bits)].get();
        else
          return n->value;
    }
  };
}

char hex_digit(Int a)
{
  Int b = a & 0xf;
  return b < 10 ? '0' + b : 'A' + (b - 10);
}

void test_encode_decode(const std::vector<Byte>& bytes, int n)
{
  std::vector<Int> freq(n, 0);
  for(auto e : bytes) freq[e]++;

  std::vector<Int> encode_tree;
  huffman::create_encode_tree(&freq[0], freq.size(), encode_tree);

  BitArray encoded;
  huffman::Encoder encode(encode_tree);
  for(auto e : bytes)
  {
    auto tmp = encode(e);
    encoded.push(tmp);
    cout << "e " << Int(e) << endl;
  }

  huffman::Decoder decode(encode_tree);
   
  BitView v = encoded.get_view();
  std::vector<Int> decoded;
  while(decoded.size() < bytes.size())
    decoded.push_back(decode(v));

  Int ndiff = 0;
  for(Int i = 0; i < bytes.size(); i++)
    if(bytes[i] != decoded[i])
      ndiff++;

  for(Int i = 0; i < freq.size(); i++)
    cout << "freq " << i << " " << freq[i] << endl;

  for(Int i = 0; i < n; i++)
  {
    cout << i << ": ";
    cout << encode.table[i].offset << " ";
    for(auto e : encode.table[i].v) cout << hex_digit(e) << hex_digit(e >> 4);
    cout << " ";
    for(auto v = encode(i); v.len > 0; ) cout << Int(read_bit(v));
    cout << endl;
  }
  
  for(Int i = 0; i < n; i++)
  {
  }

  cout << "num diff " << ndiff << endl;
  cout << "decoded size " << decoded.size() << endl;
  cout << "encoded size " << encoded.v.size() << endl;

}

int main(int argc, char** argv)
{
  if(argc != 4) return 1;

  if(strcmp(argv[1], "e") == 0)
  {
    Int n = 256;

    auto img = load_pgm(argv[2]);
    std::vector<Byte> encoded;
    runlength_encode(img.ptr.get(), img.sz0 * img.sz1, n - 1, encoded);

    test_encode_decode(encoded, n);

    uint16_t sz0 = img.sz0;
    uint16_t sz1 = img.sz1;
    uint32_t len = encoded.size();

    std::ofstream f(argv[3], std::ios_base::binary);
    f.write((const char*) &sz0, sizeof(sz0));
    f.write((const char*) &sz1, sizeof(sz1));
    f.write((const char*) &len, sizeof(len));
    f.write((const char*) &encoded[0], len);
  }
  else
  {
    uint16_t sz0;
    uint16_t sz1;
    uint32_t len;
    std::ifstream f(argv[2], std::ios_base::binary);
    f.read((char*) &sz0, sizeof(sz0)); 
    f.read((char*) &sz1, sizeof(sz1));
    f.read((char*) &len, sizeof(len));
    std::vector<Byte> encoded(len);
    f.read((char*) &encoded[0], len);

    std::vector<Byte> decoded;
    runlength_decode(&encoded[0], encoded.size(), decoded);

    if(decoded.size() != sz0 * sz1) return 1;

    auto img = Img<Byte>(sz0, sz1);
    std::copy(decoded.begin(), decoded.end(), img.ptr.get());
    save_pgm(img, argv[3]);    
  }

  return 0;
}
