import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import torch
import torch.nn as nn

def combined_loss_plot(epoch, pl_dirname, point_x, point_pi, *loss_data, save=True):
    def loss_plot(ax1, epoch_num, train_loss_list, val_loss_list, T_list, title):
        epoch_num_train = range(0, epoch_num+1, 1)
        epoch_num_val = range(0, epoch_num+1, 10)

        ax2 = ax1.twinx()
        if title != 'elbo':
            ax2.plot(epoch_num_train, T_list, label='Temperature', c='0.6')
            ax2.set_ylabel('log_Temperature')
            ax2.set_yscale('log')
        ax1.plot(epoch_num_train, train_loss_list, label='train')
        ax1.plot(epoch_num_val, val_loss_list, label='val')
        ax1.axvline(point_x, ls="--", color="m")
        ax1.axvline(point_pi, ls="--", color="c")
        ax1.grid()

        handler1, label1 = ax1.get_legend_handles_labels()
        handler2, label2 = ax2.get_legend_handles_labels()
        ax2.legend(handler1 + handler2, label1 + label2, loc='upper right', borderaxespad=0.)

        ax1.set_ylabel('log_loss')
        ax1.set_yscale('log')
        ax1.set_title(title)

    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    titles = ['elbo', 'X', 'B', 'KL', 'sp_pi', 'sp_phi']

    for ax, data, title in zip(axes.ravel(), loss_data, titles):
        train_loss, loss, T_list = data
        loss_plot(ax, epoch, train_loss, loss, T_list, title)

    if save==True:
        plt.tight_layout()
        fig.savefig(pl_dirname + f'/loss_plots.pdf')
        plt.close(fig)

def π_learn_plot(epoch_num, pi_list, dirname, save=True):
    point_x = epoch_num['train_x']
    point_pi = epoch_num['train_pi']
    epoch_num_val = range(0,epoch_num['train_sp'],10)
    fig = plt.figure()
    ax = fig.add_subplot()
    for i in range(0, len(pi_list)):
        pi_list[i] = pi_list[i]
        ax.plot(epoch_num_val,np.array(pi_list[i][:]).mean(axis=1),label=f'{i}')
    ax.axvline(point_x, ls = "--", color = "m")
    ax.axvline(point_pi, ls = "--", color = "c")
    #ax.set_title(f'lr = {lr}')
    ax.set_xlabel('epoch')
    ax.set_ylabel('probability')
    ax.grid()
    if save == True:
        fig.savefig(dirname + f'/pi_learn_plots.pdf')
        plt.close(fig)


def π_plot(pi_sort, dirname, save=True):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.bar(x=range(pi_sort.shape[1]), height=(pi_sort*100).mean(axis=0), tick_label=list(range(pi_sort.shape[1])))
    #ax.set_title(f'lr = {sp_lambda}')
    ax.set_xlabel('clone')
    ax.set_ylabel('probability')
    for i in range(1, pi_sort.shape[1]+1):
        ax.text(i-1,(pi_sort*100).mean(axis=0)[i-1]+0.005,"%.2f" % (pi_sort*100).mean(axis=0)[i-1], ha='center', fontsize="xx-small")
    if save == True:
        fig.savefig(dirname + f'/pi_plots.pdf')
        plt.close(fig)

def phi_plot(phi, phi0, phik, alpha, dirname, mutations, save=True):
    fig, ax = plt.subplots(2,3, figsize=(20,10))
    ax = ax.flatten()
    sns.heatmap(phi, ax=ax[0], square=True, yticklabels=mutations).set_title('phi_sort')
    sns.heatmap(alpha, ax=ax[1], square=True, yticklabels=mutations).set_title('alpha_common')
    sns.heatmap(phi0, ax=ax[2], square=True, yticklabels=mutations).set_title('phi_common')
    sns.heatmap(1-alpha, ax=ax[4], square=True, yticklabels=mutations).set_title('alpha_specific')
    sns.heatmap(phik, ax=ax[5], square=True, yticklabels=mutations).set_title('phi_specific')
    
    for i in range(6):
        ax[i].set_xlabel('clone')
        ax[i].set_ylabel('mutation')
    ax[3].axis('off')

    if save == True:
        fig.savefig(dirname+f'/phi_plots.pdf')
        plt.close(fig)


def embedding_plot(embedding, adata, adata_mg, dirname, cell_label, mutation_label, palette=None, s=2, save=True):
    # mutation_label(genotype)/cell_label(celltype) が None または adata.obs に無い場合は
    # そのパネルをスキップする。ATAC など genotype アノテーションが無いデータに対応。
    labels = [lab for lab in (mutation_label, cell_label) if lab is not None and lab in adata.obs.columns]
    if len(labels) == 0:
        return

    fig, ax = plt.subplots(1, len(labels), figsize=(5*len(labels), 3), squeeze=False)
    ax = ax.flatten()
    for i, lab in enumerate(labels):
        hue = adata.obs[lab].astype('category')
        palette_used = palette if palette is not None else ('tab10' if len(hue.cat.categories) < 10 else 'tab20')
        sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], s=s, hue=hue, palette=palette_used, ax=ax[i]).set_title(lab)
        ax[i].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0, fontsize='x-small')
        ax[i].axes.xaxis.set_visible(False)
        ax[i].axes.yaxis.set_visible(False)
        ax[i].set_aspect('equal', adjustable='datalim')
    if save==True:
        fig.savefig(dirname + f'/embedding_plots.pdf', bbox_inches='tight')
        plt.close(fig)


def pi_emb_plot(embedding, pi_sorted, dirname, s=2, save=True):
    num_categories = pi_sorted.shape[1]
    n_rows = int(np.ceil(num_categories / np.sqrt(num_categories)))
    n_cols = n_rows

    fig, ax = plt.subplots(n_rows, n_cols, figsize=(7*n_cols, 4*n_rows))
    if n_rows == 1 and n_cols == 1:
        ax = np.array([ax])
    else:
        ax = ax.flatten()

    cmap = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)

    for i in range(num_categories):
        # データをソート
        values = pi_sorted[:,i]
        sort_indices = np.argsort(values)
        sorted_x = embedding[sort_indices, 0]
        sorted_y = embedding[sort_indices, 1]
        sorted_values = values[sort_indices]

        scatter = ax[i].scatter(
            sorted_x,
            sorted_y,
            c=sorted_values,
            s=s,
            cmap=cmap
        )
        ax[i].set_title(f'clone{i}_probability(π)', fontsize=16)
        plt.colorbar(scatter, ax=ax[i])
        ax[i].axes.xaxis.set_visible(False)
        ax[i].axes.yaxis.set_visible(False)
        ax[i].set_aspect('equal', adjustable='datalim')

    for j in range(num_categories, len(ax)):
        ax[j].axis('off')

    plt.tight_layout()

    if save:
        fig.savefig(dirname + f'/emb_pi_plots.pdf')
        plt.close(fig)