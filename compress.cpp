#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <queue>
#include <unordered_map>

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
  
  void push(Int val, Int n) { for(Int i = 0; i < n; i++) push(val >> i); }

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

void runlength_encode(const Byte* src, Int n, Int max_len, std::vector<Int>& dst)
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

void runlength_decode(const Int* src, ptrdiff_t n, std::vector<Byte>& dst)
{
  for(Int i = 0; i < n; i += 2)
    std::fill_n(std::back_inserter(dst), src[i + 1], src[i]);
}

Int wrap_lo(Int a, Int len) { return a >= 0 ? a : a + len; }
Int wrap_hi(Int a, Int len) { return a < len ? a : a - len; }

template<typename T>
void encode_img(const Img<T>& img, std::vector<Int>& dst)
{
  dst.push_back(img.sz0);  
  dst.push_back(img.sz1);
  Img<bool> visited(img.sz0, img.sz1);
  Int prev_i0 = 0;
  Int prev_i1 = 0;
  for_img(img, i0, i1)
    if((i1 == 0 || img(i0, i1) != img(i0, i1 - 1)) && !visited(i0, i1))
    {
      Int color = img(i0, i1);
      visited(i0, i1) = true;
      dst.push_back(img.sz1 + color);
      dst.push_back(i0 - prev_i0);
      dst.push_back(wrap_lo(i1 - prev_i1, img.sz1));
      prev_i0 = i0;
      prev_i1 = i1;

      Int j1 = i1;
      for(Int j0 = i0 + 1; j0 < img.sz0; j0++)
      {
        Int prev_j1 = j1;
        if(img(j0, j1) == color)
        {
          for(; j1 > 0 && img(j0, j1 - 1) == color; j1--) { } 

          if(!visited(j0, j1))
          {
            visited(j0, j1) = true;
            dst.push_back(wrap_lo(j1 - prev_j1, img.sz1));
          }
          else break;
        }
        else
        {
          for(
            j1++;
            j1 < img.sz1 && img(j0 - 1, j1) == color && img(j0, j1) != color;
            j1++) { }

          if(j1 < img.sz1 && img(j0, j1) == color && !visited(j0, j1))
          {
            visited(j0, j1) = true;
            dst.push_back(wrap_lo(j1 - prev_j1, img.sz1));
          }
          else break;
        }
      } 
    }
}

Img<Byte> decode_img(Int* ptr, Int len)
{
  Int* end = ptr + len;
  Int sz0 = *ptr; ptr++;
  Int sz1 = *ptr; ptr++;
  Int i0 = 0;
  Int i1 = 0;
  Img<Byte> img(sz0, sz1);
  Int color = *ptr - sz1; ptr++;
  while(ptr < end)
  {
    i0 = i0 + *ptr; ptr++;
    i1 = wrap_hi(i1 + *ptr, sz1); ptr++;
    Int j0 = i0;
    Int j1 = i1; 
    img(j0, j1) = color;
    while(ptr < end)
    {
      int a = *ptr; ptr++;
      if(a >= sz1)
      {
        color = a - sz1;
        break;
      }

      j1 = wrap_hi(j1 + a, sz1);
      j0++;
      img(j0, j1) = color;
    }
  }

  return img;
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

  Int num_used_bits(Int val)
  {
    Int nbits = 0;
    for(auto a = val; a; a >>= 1) nbits++;
    return nbits;
  }

  void encode(Int* ptr, Int len, std::vector<Byte>& dst)
  {
    BitArray bits;
    Int num_values = *std::max_element(ptr, ptr + len) + 1;
    Int nbits = num_used_bits(num_values);
    bits.push(nbits, 8);

    int n;
    std::vector<Int> freq;
    std::vector<Int> table;
    {
      std::vector<std::pair<Int, Int>> tmp;
      for(Int i = 0; i < num_values; i++) tmp.emplace_back(0, i);
      for(Int i = 0; i < len; i++) tmp[ptr[i]].first++;
      std::sort(tmp.begin(), tmp.end(), std::greater<>());
      Int sum = 0; for(auto e : tmp) sum += e.first;
      Int log2_sum = num_used_bits(sum);  

      Int n = 0;
      bool n_initialized = false;
      Int n_sz = 0;

      Int remaining_sum = sum;
      Int num_bits_special = 0;
      for(Int i = 0;; i++)
      {
        int current_n = i + 1;

        Int num_bits_general =
          remaining_sum * (log2_sum - num_used_bits(remaining_sum) + nbits);

        Int tree_nbits = num_used_bits(i);
        Int sz = 
          num_bits_general +
          num_bits_special +
          (current_n - 1) * nbits +
          (2 * current_n - 1) * tree_nbits;

        if(!n_initialized || sz < n_sz)
        {
          n_sz = sz;
          n = current_n;
          n_initialized = true;
        }

        if(remaining_sum == 0 || i == tmp.size()) break;

        Int freq = tmp[i].first;
        num_bits_special += (log2_sum - num_used_bits(freq)) * freq;
        remaining_sum -= freq;
      }

      for(Int i = 0; i < n - 1; i++)
      {
        freq.push_back(tmp[i].first);
        table.push_back(tmp[i].second);
      }

      freq.push_back(0);
      for(Int i = n - 1; i < num_values; i++) freq.back()++;
    }

    std::unordered_map<Int, Int> reverse_table;
    for(Int i = 0; i < table.size(); i++) reverse_table.emplace(table[i], i);

    std::vector<Int> encode_tree;
    huffman::create_encode_tree(&freq[0], freq.size(), encode_tree);

    bits.push(n, nbits);
    for(auto e : table) bits.push(e, nbits);

    Int tree_nbits = num_used_bits(n);
    for(Int i = 0; i < encode_tree.size() - 1; i++)
      bits.push(encode_tree[i] - n, tree_nbits);

    Encoder encode(encode_tree);
    bits.push(len, 32);
    for(Int i = 0; i < len; i++)
    {
      Int val = ptr[i];
      auto it = reverse_table.find(val);
      if(it == reverse_table.end())
      {
        bits.push(encode(table.size()));
        bits.push(val, nbits);
      }
      else bits.push(encode(it->second));
    }

    dst = std::move(bits.v);
  }
  
  void decode(Byte* ptr, Int len, std::vector<Int>& dst)
  {
    BitView bits = { ptr, 0, len * bytesize };
    Int nbits = read_bits(bits, 8); 
    Int n = read_bits(bits, nbits);
    std::vector<Int> table;
    for(Int i = 0; i < n - 1; i++) table.push_back(read_bits(bits, nbits));
    std::vector<Int> encode_tree;
    Int tree_nbits = num_used_bits(n);
    for(Int i = 0; i < 2 * n - 2; i++)
      encode_tree.push_back(n + read_bits(bits, tree_nbits));

    encode_tree.push_back(-1);

    Decoder decode(encode_tree);
    Int dst_len = read_bits(bits, 32);
    dst.clear();
    for(Int i = 0; i < dst_len; i++)
    {
      Int a = decode(bits);
      if(a == table.size())
        dst.push_back(read_bits(bits, nbits));
      else
        dst.push_back(table[a]);
    }
  }
}

char hex_digit(Int a)
{
  Int b = a & 0xf;
  return b < 10 ? '0' + b : 'A' + (b - 10);
}

void test_encode_decode1(const std::vector<Byte>& bytes, int n)
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
    for(auto v = encode(i); v.len > 0; ) cout << Int(read_bit(v));
    cout << endl;
  }
  
  cout << "num diff " << ndiff << endl;
  cout << "decoded size " << decoded.size() << endl;
  cout << "encoded size " << encoded.v.size() << endl;
}

int main(int argc, char** argv)
{
  typedef std::chrono::duration<double> Dur;

  if(argc != 4) return 1;

  if(strcmp(argv[1], "e") == 0)
  {
    Int n = 1 << 8;

    auto img = load_pgm(argv[2]);
    std::vector<Int> intermediate;

    auto start = std::chrono::system_clock::now();

#if 0
    runlength_encode(img.ptr.get(), img.sz0 * img.sz1, n - 1, intermediate);
#else
    encode_img(img, intermediate);
#endif

    auto end_ei = std::chrono::system_clock::now();

    std::vector<Byte> huff_encoded;
    //for(Int i = 0; i < 1000; i++)
    huffman::encode(&intermediate[0], intermediate.size(), huff_encoded);

    auto end_huff = std::chrono::system_clock::now();

    cout <<
      Dur(end_ei - start).count() << " " <<
      Dur(end_huff - end_ei).count() << endl;

    uint16_t sz0 = img.sz0;
    uint16_t sz1 = img.sz1;
    uint32_t len = huff_encoded.size();
    std::ofstream f(argv[3], std::ios_base::binary);
    f.write((const char*) &sz0, sizeof(sz0));
    f.write((const char*) &sz1, sizeof(sz1));
    f.write((const char*) &len, sizeof(len));
    f.write((const char*) &huff_encoded[0], len);
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
    std::vector<Byte> huff_encoded(len);
    f.read((char*) &huff_encoded[0], len);

    std::vector<Int> intermediate(len);
    //for(Int i = 0; i < 1000; i++)
    huffman::decode(&huff_encoded[0], huff_encoded.size(), intermediate);

    std::vector<Byte> decoded;

#if 0
    runlength_decode(&intermediate[0], intermediate.size(), decoded);
    if(decoded.size() != sz0 * sz1) return 1;
    auto img = Img<Byte>(sz0, sz1);
    std::copy(decoded.begin(), decoded.end(), img.ptr.get());
#else
    auto img = decode_img(&intermediate[0], intermediate.size());
#endif

    save_pgm(img, argv[3]);    
  }

  return 0;
}
