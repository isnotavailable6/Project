### **project4.cpp 补充说明文档**
1. **SM3 SIMD加速**
   - `ROTL_SIMD(__m256i x, int n)`：利用AVX2指令集实现32位整数的并行循环左移，支持8个数据同时处理，用于SM3压缩函数的位运算加速。
   - `P0_SIMD(__m256i X)`/`P1_SIMD(__m256i X)`：SIMD版本的SM3置换函数`P0`和`P1`，通过并行移位与异或提升运算效率。
   - `FF_SIMD(__m256i x, __m256i y, __m256i z, int j)`/`GG_SIMD(...)`：SIMD版本的布尔函数`FF`和`GG`，根据轮数`j`选择不同逻辑运算，支持8组数据并行处理。
   - `CF_SIMD(__m256i V[8], uint8_t* block)`：SM3压缩函数的SIMD实现，一次处理8个64字节消息块，通过AVX2指令并行执行64轮迭代。
   - `SM3_encrypt_SIMD(const std::vector<std::vector<uint8_t>>& messages, std::vector<std::vector<uint8_t>>& hash_outputs)`：批量处理8个64字节消息的SM3哈希函数，基于`CF_SIMD`实现并行计算。

2. **Merkle树扩展**
   - `rfc_leaf_hash(const vector<uint8_t>& data)`：按规范计算叶子节点哈希，在原始数据前添加`0x00`前缀后进行SM3哈希。
   - `rfc_node_hash(const vector<uint8_t>& left, const vector<uint8_t>& right)`：按规范计算内部节点哈希，在左右子节点哈希前添加`0x01`前缀后进行SM3哈希。
   - `MerkleTree::get_inclusion_proof(int index)`：生成指定索引叶子节点的包含性证明，收集路径上所有兄弟节点的哈希值。
   - `verify_inclusion_proof(...)`：验证叶子节点的包含性证明，通过迭代合并当前哈希与兄弟节点哈希，最终对比根哈希验证有效性。

3. **其他辅助函数**
   - `from_32_to_8(uint32_t mess, uint8_t* b)`：将32位整数转换为4字节数组（大端序）。
   - `from_8_to_32(const uint8_t* b)`：将4字节数组（大端序）转换为32位整数。
   - `SM3_length_extension_attack(...)`：实现SM3哈希长度扩展攻击，基于已知哈希值和原始数据长度，构造新的哈希值对应扩展后的数据。
