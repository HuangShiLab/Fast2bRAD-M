use anyhow::{Context, Result, anyhow};
use csv::{ReaderBuilder, WriterBuilder};
use ndarray::Array1; 
use ort::value::Value;
// 修复：针对 ort 2.x 准确的导入路径
use ort::session::{builder::GraphOptimizationLevel, Session}; 
use std::path::Path;

/// 计算特征 (保持不变)
fn calculate_features(
    s_tags: f64, 
    s_reads: f64, 
    t_tags: f64, 
    total_reads_sum: f64
) -> [f32; 4] {
    let safe_s_tags = s_tags.max(1.0);
    let safe_t_tags = t_tags.max(1.0);
    let safe_s_reads = s_reads.max(1.0);
    let safe_total = total_reads_sum.max(1.0);

    let f1 = (safe_s_tags / safe_t_tags).ln();
    let g_score = (safe_s_tags * safe_s_reads).sqrt();
    let f2 = g_score.ln();
    let f3 = (safe_s_reads / safe_s_tags).ln();
    let avg_depth = safe_s_reads / safe_s_tags;
    let theoretical_reads = avg_depth * safe_t_tags;
    let f4 = (theoretical_reads / safe_total).ln();

    [f1 as f32, f2 as f32, f3 as f32, f4 as f32]
}

pub fn run_prediction(
    input_file: &Path, 
    output_file: &Path, 
    model_path: &Path
) -> Result<()> {
    // 1. 初始化 Session (必须为 mut 才能调用 run)
    // 【关键修复】使用 .map_err 将非线程安全的 ort 错误转为 anyhow 兼容的字符串错误
    let mut session = Session::builder()
        .map_err(|e| anyhow!("SessionBuilder 错误: {:?}", e))?
        .with_optimization_level(GraphOptimizationLevel::Level3)
        .map_err(|e| anyhow!("设置优化等级失败: {:?}", e))?
        .with_intra_threads(1)
        .map_err(|e| anyhow!("设置线程数失败: {:?}", e))?
        .commit_from_file(model_path)
        .map_err(|e| anyhow!("模型加载失败: {:?}", e))?;

    // 2. 读取输入文件
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input_file)
        .with_context(|| format!("无法读取输入文件: {:?}", input_file))?;
    
    let headers = reader.headers()?.clone();
    let idx_s_tags = headers.iter().position(|h| h == "Sequenced_Tag_Num").context("列丢失")?;
    let idx_s_reads = headers.iter().position(|h| h == "Sequenced_Reads_Num").context("列丢失")?;
    let idx_t_tags = headers.iter().position(|h| h == "Theoretical_Tag_Num").context("列丢失")?;

    let mut records = Vec::new();
    let mut total_reads_sum = 0.0;
    for result in reader.records() {
        let record = result?;
        total_reads_sum += record[idx_s_reads].parse::<f64>().unwrap_or(0.0);
        records.push(record);
    }
    if records.is_empty() { return Ok(()); }

    // 3. 构建特征
    let mut features_flat = Vec::with_capacity(records.len() * 4);
    for rec in &records {
        let feats = calculate_features(
            rec[idx_s_tags].parse().unwrap_or(0.0),
            rec[idx_s_reads].parse().unwrap_or(0.0),
            rec[idx_t_tags].parse().unwrap_or(1.0),
            total_reads_sum
        );
        features_flat.extend_from_slice(&feats);
    }

    let n_samples = records.len();
    
    // 【关键修复】使用元组格式 (维度, 扁平数据) 创建 Value，避免 ndarray 版本不匹配问题
    let input_value = Value::from_array(([n_samples, 4], features_flat.into_boxed_slice()))
        .map_err(|e| anyhow!("创建输入张量失败: {:?}", e))?;

    // 4. 执行推理 (inputs! 宏在 rc.12 中不返回 Result，不要加 ?)
    let outputs = session.run(ort::inputs!["float_input" => input_value])
        .map_err(|e| anyhow!("推理运行失败: {:?}", e))?;
    
    // 5. 提取结果 (处理 try_extract_tensor 返回的元组)
    let (_shape, label_data) = outputs[0].try_extract_tensor::<i64>()
        .map_err(|e| anyhow!("结果提取失败: {:?}", e))?;
    let labels_view = Array1::from_iter(label_data.iter().copied());

    // 6. 写入结果
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(output_file)?;
    let mut new_headers = headers.clone();
    new_headers.push_field("Prediction");
    writer.write_record(&new_headers)?;

    for (i, rec) in records.iter().enumerate() {
        let mut new_rec = rec.clone();
        new_rec.push_field(&labels_view[i].to_string());
        writer.write_record(&new_rec)?;
    }

    Ok(())
}