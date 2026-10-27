import torch
import torch.distributions as dist
import torch.nn as nn
#from torch import Functional as F
import numpy as np
from funcs import calc_kld, calc_nb_loss, calc_binom_cat_loss, calc_binom_cat_loss_vectorized, calc_sparce_pi, calc_sparce_phi


class LinearReLU(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(LinearReLU, self).__init__() # super().__init__() 継承元のコンストラクタ（__init()__）をオーバーライド
        self.f = nn.Sequential(
            #nn.Dropout(p=0.1),
            nn.Linear(input_dim, output_dim), # y = x A^T + b の x にデータが入る
            nn.LayerNorm(output_dim, elementwise_affine=False)) # レイヤー正則化
            # 活性化関数は置かない (線形層 + LayerNorm、非線形性は LayerNorm だけ)。
            # 旧版の nn.LeakyReLU(True) は第1引数が negative_slope のため傾き 1 の恒等関数で、論文の fit はすべてこの構成。
            # nn.LeakyReLU() (傾き 0.01) にすると CLL で段階 1 の早い時点から過学習し、z が LSI クラスタの構造を失った
            # (2026-10-10 の診断: 活性化なしは 600 エポックまで検証損失が下がり続け、LeakyReLU は 90 エポックで最良)

    def forward(self, x): # 順伝播？上のself.fを実行して出力を格納？
        h = self.f(x)
        return(h)


class SeqNN(nn.Module):
    def __init__(self, num_steps, dim):
        super(SeqNN, self).__init__()
        modules = [
            LinearReLU(dim, dim)
            for _ in range(num_steps) # _ ってなんだ？ → 使っていない変数、ただ繰り返すための変数 LinearReLUをnum_stepsだけ繰り返すってこと
        ]
        self.f = nn.Sequential(*modules)

    def forward(self, pre_h):
        post_h = self.f(pre_h)
        return(post_h)


class Encoder(nn.Module):
    def __init__(self, num_h_layers, x_dim, h_dim, z_dim):
        super(Encoder, self).__init__()
        self.x2h = LinearReLU(x_dim, h_dim) # xからzへのNN
        self.seq_nn = SeqNN(num_h_layers - 1, h_dim) #
        #self.dropout = nn.Dropout(p=0.1)
        self.h2mu = nn.Linear(h_dim, z_dim)
        self.h2logvar = nn.Linear(h_dim, z_dim)

    def forward(self, x):
        pre_h = self.x2h(x)
        post_h = self.seq_nn(pre_h)
        #post_h = self.dropout(post_h)
        mu = self.h2mu(post_h + 1e-7) # add 1e-7 to avoid nan
        logvar = self.h2logvar(post_h)
        #lib_size = torch.exp(torch.sum(x))
        return(mu, logvar)


class Decoder_1(nn.Module):
    def __init__(self, num_h_layers, z_dim, h_dim, x_dim):
        super(Decoder_1, self).__init__()
        self.z2h_ld = LinearReLU(z_dim, h_dim) #
        self.seq_nn_ld = SeqNN(num_h_layers - 1, h_dim) # hidden layer が 2 なら step は 2-1=1
        self.h2ld = nn.Linear(h_dim, x_dim) # 入力次元がh_dim, 出力次元がx_dim つまりデコード lambda の出力
        self.softplus = nn.Softplus()
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, z):
        pre_h_ld = self.z2h_ld(z)
        post_h_ld = self.seq_nn_ld(pre_h_ld)
        ld = self.h2ld(post_h_ld)
        # norm_mat (= 細胞の総数 x 特徴の総数 / 全体の総数、独立モデルの期待カウント) に掛ける正の倍率なので softplus。
        # 旧版の softmax だと特徴方向の和が 1 に縛られ、NB の平均が観測の約 1/x_dim になっていた (2026-10-09 修正)
        correct_ld = self.softplus(ld)
        #px_rate = correct_ld
        
        return(correct_ld)
    
class Decoder_2(nn.Module):
    def __init__(self, num_h_layers, z_dim, h_dim, k):
        super(Decoder_2, self).__init__()
        self.z2h_pi = LinearReLU(z_dim, h_dim) # 追加した
        self.seq_nn_pi = SeqNN(num_h_layers - 1, h_dim) # 追加した
        self.h2pi = nn.Linear(h_dim, k) # 追加した クラスタ数kどうする？
        self.softmax = nn.Softmax(dim=1) # 追加した
        self.softplus = nn.Softplus()
        self.ReLU = nn.ReLU()

    def forward(self, z):
        pre_h_pi = self.z2h_pi(z)
        post_h_pi = self.seq_nn_pi(pre_h_pi)
        pi = self.h2pi(post_h_pi) # 追加した
        pi_z = self.softmax(pi)

        return(pi_z)

class scVAE(nn.Module):
    def __init__(
            self,
            x_dim, xz_dim,# 入力と出力の次元
            enc_z_h_dim, dec_z_h_dim,# エンコーダーとデコーダーのzの次元
            num_enc_z_layers,# エンコーダーの層の数
            num_dec_z_layers, k, mutations, n_samples, use_phi_common=False, **kwargs):# デコーダーの層の数と、可変長引数（？）kを追加した
        super(scVAE, self).__init__()
        self.enc_z = Encoder(num_enc_z_layers, x_dim+n_samples, enc_z_h_dim, xz_dim)
        self.dec_z2x = Decoder_1(num_dec_z_layers, xz_dim+n_samples, dec_z_h_dim, x_dim)
        self.dec_z2pi = Decoder_2(num_dec_z_layers, xz_dim+n_samples, dec_z_h_dim, k)
        self.softplus = nn.Softplus()
        self.logtheta_x = nn.Parameter(torch.Tensor(x_dim))
        self.logphi_k = nn.Parameter(torch.Tensor(k,len(mutations)))
        self.use_phi_common = use_phi_common
        if use_phi_common:
            # バックグラウンド変異や共通変異が存在する場合に有効化
            self.logphi_0 = nn.Parameter(torch.Tensor(1,len(mutations)))
            self.alpha = nn.Parameter(torch.Tensor(k,1))
        self.sigmoid = nn.Sigmoid()
        self.reset_parameters()
        self.k = k
        self.pi_list = [[] for i in range(k)]
        #self.pi_list = torch.zeros(k, 0)

    def reset_parameters(self, seed=None):
        if seed is not None:
            torch.manual_seed(seed)
        nn.init.normal_(self.logtheta_x)
        nn.init.normal_(self.logphi_k)
        if self.use_phi_common:
            nn.init.normal_(self.logphi_0)
            nn.init.normal_(self.alpha)

    def encode_qz(self, x, sample_info):
        # q(z|x) の平均と標準偏差。forward と、encoder 凍結後のキャッシュ (exp.VaeExperiment._build_qz_cache) が共用する
        if sample_info.dim() > 2:
            sample_info = sample_info.squeeze(-1)
        x = torch.cat([x, sample_info], axis=1)
        qz_mu, qz_logvar = self.enc_z(x) # 引数一個だけでいいの？ self.enc_z = Encoder(num_enc_z_layers, x_dim, enc_z_h_dim, xz_dim)
        qz_logvar = self.softplus(qz_logvar) # これどうなんだ？
        return(qz_mu, qz_logvar)

    def calc_phi(self):
        phi_k = self.sigmoid(self.logphi_k)

        if self.use_phi_common:
            alpha = self.sigmoid(self.alpha)
            phi_0 = self.sigmoid(self.logphi_0.expand(self.logphi_k.size()))
            phi = alpha * phi_0 + (1 - alpha) * phi_k
        else:
            phi = phi_k
            phi_0 = torch.zeros_like(phi_k)
            alpha = torch.zeros(self.k, 1)
        return(phi, phi_0, phi_k, alpha)

    def forward(self, x, sample_info):
        # encode z
        if sample_info.dim() > 2:
            sample_info = sample_info.squeeze(-1)
        #print('sample_info type 1:', type(sample_info))
        #print('sample_info sum_encz', sample_info.sum(dim=0))
        qz_mu, qz_logvar = self.encode_qz(x, sample_info)
        qz = dist.Normal(qz_mu, qz_logvar) # 平均qz_mu, 偏差（分散？）はqz_logvarをソフトプラスに通した、正規分布で出力されるqz
        z = qz.rsample() # z ~ q(z|x)
        # decode z
        #print('sample_info type 2:', type(sample_info))
        #print('sample_info sum_decz', sample_info.sum(dim=0))
        z = torch.cat([z, sample_info], axis=1)
        xld = self.dec_z2x(z) # xld は、サンプリングされたzを使ってデコードされたp(x|z)
        pi_z = self.dec_z2pi(z) # xld は、サンプリングされたzを使ってデコードされたp(x|z)
        theta_x = self.softplus(self.logtheta_x)
        phi, phi_0, phi_k, alpha = self.calc_phi()

        return(z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha)

    def elbo_loss(self, x, xnorm_mat, d, b, sample_info, params):
        z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha = self(x, sample_info) # ここで、z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_kを求める
        # 旧版はここで毎バッチ pi を self.pi_list に貯めていた (π_learn_plot 用、現在は未使用)。
        # クローン数ぶんの GPU→CPU 同期と、学習の長さに比例するメモリ増加 (ovarian で 20 GB 超) になるので廃止
        elbo_loss = 0 # kldとXのみのloss train_x
        kld, kld_T = calc_kld(qz, params['kld_temp'].temp) # kldの計算
        elbo_loss += kld_T.mean()
        # reconst loss of xS
        #X_recon_loss = calc_poisson_loss(xld, xnorm_mat, x).mean()
        
        theta_x = self.softplus(self.logtheta_x) + 1e-6  # 小さな値を追加
        theta_x = torch.clamp(theta_x, min=1e-5, max=1e5)  # 安全範囲にクランプ
        X_recon_loss = calc_nb_loss(xld, xnorm_mat, theta_x, x).mean()
        elbo_loss += X_recon_loss

        B_recon_loss, B_recon_loss_T = calc_binom_cat_loss_vectorized(pi_z, d, b, phi, params['binom_temp'].temp) # B の reconstraction loss 記録

        elbo_loss_B = elbo_loss + B_recon_loss_T.mean()
        sp_loss_pi, sp_pi_T = calc_sparce_pi(pi_z, params['lmd_pi'], params['sp_temp'].temp)
        sp_loss_phi, sp_phi_T = calc_sparce_phi(phi, params['lmd_phi'], params['sp_temp'].temp)
        elbo_loss_sp = elbo_loss_B + sp_pi_T + sp_phi_T
        #elbo_loss_sp = elbo_loss_B + sp_pi_T

        return(elbo_loss, elbo_loss_B, elbo_loss_sp, kld.mean(), X_recon_loss, B_recon_loss.mean(), sp_loss_pi, sp_loss_phi)

    def elbo_loss_from_qz(self, qz_mu, qz_scale, d, b, sample_info, params):
        # encoder を凍結した train_pi / train_sp 用。x の代わりに q(z|x) の平均・標準偏差を受け取る。
        # このとき x の再構成項は学習中のパラメータ (dec_z2pi, phi) に勾配を持たないので計算しない。
        # z の引き方は forward と同じ (rsample 1 回) なので、乱数の消費も elbo_loss と一致する。
        # 戻り値は elbo_loss と同じ並び。x の再構成項は計算していないので NaN を入れ、合計には含めない
        # qz_scale=None (外部の埋め込みなど) なら z を平均に固定し、KL は定義されないので 0 とする (学習中のパラメータに勾配を持たない)
        if sample_info.dim() > 2:
            sample_info = sample_info.squeeze(-1)
        if qz_scale is None:
            z = qz_mu
            kld = torch.zeros(qz_mu.shape[0], device=qz_mu.device); kld_T = kld
        else:
            qz = dist.Normal(qz_mu, qz_scale)
            z = qz.rsample() # z ~ q(z|x)
            kld, kld_T = calc_kld(qz, params['kld_temp'].temp)
        z = torch.cat([z, sample_info], axis=1)
        pi_z = self.dec_z2pi(z)
        phi, phi_0, phi_k, alpha = self.calc_phi()

        elbo_loss = kld_T.mean()
        X_recon_loss = torch.full((), float('nan'), device=qz_mu.device)

        B_recon_loss, B_recon_loss_T = calc_binom_cat_loss_vectorized(pi_z, d, b, phi, params['binom_temp'].temp)

        elbo_loss_B = elbo_loss + B_recon_loss_T.mean()
        sp_loss_pi, sp_pi_T = calc_sparce_pi(pi_z, params['lmd_pi'], params['sp_temp'].temp)
        sp_loss_phi, sp_phi_T = calc_sparce_phi(phi, params['lmd_phi'], params['sp_temp'].temp)
        elbo_loss_sp = elbo_loss_B + sp_pi_T + sp_phi_T

        return(elbo_loss, elbo_loss_B, elbo_loss_sp, kld.mean(), X_recon_loss, B_recon_loss.mean(), sp_loss_pi, sp_loss_phi)