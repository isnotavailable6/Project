### **project1.cpp 补充说明文档**
1. **密钥派生函数**
   - `Kdf(z: bytes, klen: int)`：基于SM3哈希的密钥派生函数，通过迭代哈希`z + 计数器`生成指定长度（`klen`）的密钥材料，计数器从1开始递增，最终截取前`klen`字节作为结果。

2. **AES-NI加速相关**
   - `L_NI(uint32_t B)`：利用AES-NI指令集实现SM4算法中的线性变换`L`，通过SIMD指令并行处理32位数据的循环移位与异或操作，提升运算效率。
   - `T_NI(uint32_t x)`：结合S盒变换（`S(x)`）与AES-NI加速的线性变换`L_NI`，构成SM4轮函数中的非线性变换`T`的硬件加速版本。
   - `EncryptBlock_AES(const uint8_t input[16], uint8_t output[16], const uint32_t rk[32])`：使用AES-NI指令集加速的SM4分组加密函数，通过`T_NI`函数提升轮迭代效率。

3. **GCM模式辅助函数**
   - `GF128_rs(uint8_t V[16])`：伽罗瓦域`GF(2^128)`中的右移操作，配合多项式-1多项式实现域内乘法，用于GMAC认证标签计算。
   - `GF128_table_generate(const uint8_t H[16])`：预计算伽罗瓦域乘法表（`H_table_high`和`H_table_low`），加速GMAC中的乘法运算。
   - `GF128_Mul(uint8_t X[16], uint8_t out[16])`：利用预计算表实现`GF(2^128)`中的乘法，用于更新GMAC的中间结果。
   - `H_generate(uint8_t key[16], uint8_t H[16])`：生成GCM模式中的哈希密钥`H`，通过加密全零块得到。

4. **测试函数扩展**
   - `AES_NI_check()`：测试AES-NI加速版本的SM4加密性能，对比普通实现的耗时，验证硬件加速效果。
