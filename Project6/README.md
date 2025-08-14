### **project6.py 补充说明文档**
#### 协议流程扩展
1. **参与方1交互**
   - `Round1(self) -> List[int]`：将本地数据集的哈希值（经`k1`加密）打乱后发送给参与方2。
   - `Round3(self, pairs, z_set: set[int], pk: paillier.PaillierPublicKey)`：接收参与方2的第二轮数据，筛选出交集元素对应的加密值并求和，返回同态加密的和。

2. **参与方2交互**
   - `Round2(self, H: List[int])`：对接收到的哈希值进行`k2`加密并打乱，同时计算本地数据集的哈希值（经`k2`加密）和键值对的同态加密，返回处理后的数据。
   - `decrypt_sum(self, ciphertext)`：使用Paillier私钥解密密文，得到交集元素的求和结果。
