fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharey='row')

for n_neighbors in (3, 5, 7, 9):
  df_n = res_df[res_df.n_neighbors == n_neighbors]

  axs[0].plot(df_n.nsamples, df_n.train_acc, label=f"{n_neighbors}")
  axs[0].set_title("Train")

  axs[0].legend(loc='best')


  axs[1].plot(df_n.nsamples, df_n.test_acc, label=f"{n_neighbors}")
  axs[1].set_title("Test")
  axs[1].legend(loc='best')

  axs[0].set_xlabel("N samples")
  axs[1].set_xlabel("N samples")
  axs[0].set_ylabel("Accuracy")
  axs[1].set_ylabel("Accuracy")

  axs[0].grid(True)
  axs[1].grid(True)