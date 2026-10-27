import torch
import torch.distributions as dist
import torch.nn as nn
import numpy as np


# loss function
def calc_kld(qz, T):
    kld = -0.5 * (1 + qz.scale.pow(2).log() - qz.loc.pow(2) - qz.scale.pow(2)) # qz.loc = 平均, qz.scale = 標準偏差 解読完了。詳しくは21/11/11 のノートへ
    kld_sum = kld.sum(-1)
    kld_sum_T = kld_sum / T
    return(kld_sum, kld_sum_T)


def calc_nb_loss(ld, norm_mat, theta, obs, eps=1e-8):
    ld = norm_mat * ld
    # p = ld / (ld + theta) を確率で渡すと、ld が theta の約 1.7e7 倍を超えたとき float32 で p が 1 に丸まり、
    # NegativeBinomial の引数検査で例外になる。log(p / (1 - p)) = log(ld) - log(theta) なので logits で渡す
    # (log_prob は内部で logits を使うので式は同じ。eps は ld = 0 のとき log が -inf になるのを防ぐ。2026-10-10)
    logits = torch.log(ld + eps) - torch.log(theta + eps)
    p_z = dist.NegativeBinomial(theta, logits=logits)
    l = - p_z.log_prob(obs)
    l_sum = l.sum(-1)
    return(l_sum)


def calc_poisson_loss(ld, norm_mat, obs): # Poissonよりも、他の分布の方が良さそう
    p_z = dist.Poisson(ld * norm_mat) # なんでパラメータが ld * norm_mat なんだ？→norm_matは各細胞におけるexpressionの平均→情報をあらかじめ与えて学習効率を上げる
    l = - p_z.log_prob(obs)
    l_sum = l.sum(-1)
    return(l_sum)

def calc_binom_cat_loss(pi_z, d, b, phi, T): # 追加した部分
    l_pi = pi_z.log() # 変更 log pi と同じ
    l_b_sum = torch.zeros(pi_z.size(axis=-1),len(d), device=pi_z.device) # 変更

    '''
    if torch.isnan(phi).any()==True: # alphaがnanになるエラーに対してのとりあえずの対処
        phi = torch.zeros_like(phi)
    '''
    
    phi_logit  = torch.special.logit(phi)

    #scale_factor = 1.0 / torch.sqrt(torch.tensor(len(d), dtype=torch.float32))
    for k in range(pi_z.size(axis=-1)):
      p_b = dist.Binomial(total_count=d, logits=phi_logit[k]) # 確率phiで変異が入るbinomial マッピングされるリード数がd
      l_b = p_b.log_prob(b) 
      #l_b_sum[k] = scale_factor * l_b.sum(-1)
      l_b_sum[k] = l_b.sum(-1)
    
    #if l_pi.is_cuda ==True:
    #    l_b_sum = l_b_sum.to(device='cuda')

    #l_pi_b_sum = l_pi + l_b_sum.transpose(0,-1)/T
    l_pi_b = l_pi + l_b_sum.transpose(0,-1)

    # logsumexp
    l_pi_b_sum = - torch.logsumexp(l_pi_b, -1)
    l_pi_b_sum_T = - torch.logsumexp(l_pi_b/T, -1)
    return(l_pi_b_sum, l_pi_b_sum_T)

def calc_binom_cat_loss_vectorized(pi_z, d, b, phi, T):
    """
    ベクトル化された二項分布カテゴリカル損失計算関数
    
    Parameters:
    -----------
    pi_z : torch.Tensor[batch_size, num_clones]
        クローン確率
    d : torch.Tensor[batch_size, num_mutations]
        総リード数
    b : torch.Tensor[batch_size, num_mutations]
        変異リード数（b ≤ d）
    phi : torch.Tensor[num_clones, num_mutations]
        変異確率
    T : float
        温度パラメータ
    """
    # クローン確率の対数 [batch_size, num_clones]
    pi_z = pi_z.clamp(min=1e-6, max=1-1e-6)  # 数値安定性のためにクランプ
    l_pi = pi_z.log()
    
    # 数値安定性のためのphiのクランプとロジット変換
    phi_clamped = torch.clamp(phi, 1e-6, 1-1e-6)
    phi_logit = torch.special.logit(phi_clamped)
    
    # 次元拡張による3Dテンソルの準備
    d_expanded = d.unsqueeze(1)        # [batch_size, 1, num_mutations]
    b_expanded = b.unsqueeze(1)        # [batch_size, 1, num_mutations]
    phi_logit_expanded = phi_logit.unsqueeze(0)  # [1, num_clones, num_mutations]
    
    # 一度にすべての組み合わせの対数確率を計算
    # ブロードキャストにより [batch_size, num_clones, num_mutations] の形状になる
    log_probs = torch.distributions.binomial.Binomial(
        total_count=d_expanded,
        logits=phi_logit_expanded
    ).log_prob(b_expanded)
    
    # 変異にわたる合計 -> [batch_size, num_clones]
    l_b_sum = log_probs.sum(dim=2)
    
    # クローン確率の対数と合計 -> [batch_size, num_clones] 
    l_pi_b = l_pi + l_b_sum
    
    # logsumexpによるソフトマックス計算
    l_pi_b_sum = -torch.logsumexp(l_pi_b, dim=1)      # [batch_size]
    l_pi_b_sum_T = -torch.logsumexp(l_pi_b/T, dim=1)  # [batch_size]
    
    return l_pi_b_sum, l_pi_b_sum_T

def calc_sparce_pi(pi, lmd, T):
    scaled_lmd = lmd / torch.sqrt(torch.tensor(pi.size(dim=1), dtype=torch.float32))
    sum_cell = scaled_lmd * (pi.sum(dim=0)**0.5)
    sparce_term = sum_cell.sum()
    sparce_term_T = (sum_cell/T).sum()
    return(sparce_term, sparce_term_T)

def calc_sparce_phi(phi, lmd, T):
    scaled_lmd = lmd / torch.sqrt(torch.tensor(phi.size(dim=1), dtype=torch.float32))
    sparce_term = scaled_lmd * phi.sum()
    sparce_term_T = sparce_term / T
    return(sparce_term, sparce_term_T)