### **project2.py 补充说明文档**
1. **水印提取**
   - `extract_watermark(watermarked_path, original_path, watermark_shape, alpha=0.05)`：从带水印图像中提取水印。流程为：
     - 对原图和带水印图像分别进行小波分解，获取低频分量（LL）和高频分量（LH, HL, HH）；
     - 计算带水印图像与原图的LH分量差值，除以嵌入权重`alpha`得到提取的水印；
     - 调整水印尺寸至原始水印大小并返回。

2. **对比度调整**
   - `contrast(image_path)`：通过`cv2.convertScaleAbs`调整图像对比度，`alpha=1.5`增强对比度，`beta=0`保持亮度不变，结果保存为`contrast_output.png`。
