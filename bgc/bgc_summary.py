from bgc import BGC, cline_plot
import pandas as pd
import numpy as np
import os
from scipy.stats import norm

name = "clust-90-missing-0.9-exclude-1-strict" 
chains = 5 
paths = [os.path.join(name, "{}_{}.hdf5".format(name, i)) for i in range(1, chains + 1)]
rand_loci = np.sort(np.random.randint(502, size=20))

bgc = BGC(paths, 100, 0.9)

# LnL = bgc.posterior("LnL")
alpha = bgc.posterior("alpha")
tau_alpha = bgc.posterior("tau-alpha")
beta = bgc.posterior("beta")
tau_beta = bgc.posterior("tau-beta")
# hi = bgc.posterior("hi")
# gamma_quantile = bgc.posterior("gamma-quantile")
# zeta_quantile = bgc.posterior("zeta-quantile")

def plot_path(param):
    return os.path.join(name, "{}-{}.html".format(name, param))

# LnL.histogram(plot_path("LnL")) 
# alpha.histogram(plot_path("alpha"), rand_loci) 
# tau_alpha.histogram(plot_path("tau_alpha")) 
# beta.histogram(plot_path("beta"), rand_loci) 
# tau_beta.histogram(plot_path("tau_beta")) 
# hi.histogram(plot_path("hi")) 
# gamma_quantile.histogram(plot_path("gamma_quantile"), rand_loci) 
# zeta_quantile.histogram(plot_path("zeta_quantile"), rand_loci) 

# LnL.trace(plot_path("LnL")) 
# alpha.trace(plot_path("alpha"), rand_loci) 
# tau_alpha.trace(plot_path("tau_alpha")) 
# beta.trace(plot_path("beta"), rand_loci) 
# tau_beta.trace(plot_path("tau_beta")) 
# hi.trace(plot_path("hi")) 
# gamma_quantile.trace(plot_path("gamma_quantile"), rand_loci) 
# zeta_quantile.trace(plot_path("zeta_quantile"), rand_loci) 

# all_sum = pd.concat(
#     [
#         # LnL.summary,
#         alpha.summary,
#         # tau_alpha.summary,
#         beta.summary,
#         # tau_beta.summary,
#         # hi.summary,
#         # gamma_quantile.summary,
#         # zeta.summary
#     ])   
# all_sum.to_csv(os.path.join(name, "parameter-summaries.csv"))

# alpha_outliers1 = ~alpha.covered(0)
# print(f"Number of Alpha outliers: {alpha_outliers1.sum()}")

# beta_outliers1 = ~beta.covered(0)
# print(f"Number of Beta outliers: {beta_outliers1.sum()}")

# alpha_outliers = gamma_quantile.outlier(alpha.summary["median"]) 
# beta_outliers = zeta_quantile.outlier(beta.summary["median"]) 

# alpha_excess_loci = alpha.summary[alpha_excess].index.values
# alpha.histogram(plot_path("alpha-trace-excess"), alpha_excess_loci)

# cline_plot(alpha.summary["median"], beta.summary["median"], alpha_outliers)
# cline_plot(alpha.summary["median"], beta.summary["median"], beta_outliers)

from math import sqrt

def random_effects_outliers(locus, tau, interval, central_tendancy="median"):
    ## locus: posterior object for alpha or beta parameter
    ## tau: posterior object for alpha tau or alpha beta estimates
    ## interval: interval of normal probability density function 
    ## central_tendancy: "median" or "mean"
    # TODO: Not sure if this is the correct way to identify outliers
    lower_quantile = round((1.0 - interval) / 2, 3)
    upper_quantile = 1 - lower_quantile
    tau_central = tau.summary.iloc[0][central_tendancy]
    lower = norm.ppf(lower_quantile, loc=0, scale=sqrt(tau_central))
    upper = norm.ppf(upper_quantile, loc=0, scale=sqrt(tau_central))
    locus_central = locus.summary[central_tendancy]
    print(lower, upper)
    print(locus_central)
    return (locus_central < lower) | (locus_central > upper)


out = random_effects_outliers(alpha, tau_alpha, 0.8)
# print(out)
print(sum(out))

out = random_effects_outliers(beta, tau_beta, 0.8)
# print(out)
print(sum(out))














# # Linear regression
# from sklearn.linear_model import LinearRegression

# fst = pd.read_csv(os.path.join(name, name + ".fst.csv"), index_col=0)

# model = LinearRegression()
# x = fst["WEIR_AND_COCKERHAM_FST"].to_numpy()
# y = beta.summary["median"].to_numpy()

# nan = np.isnan(x)

# x = np.delete(x, nan)
# y = np.delete(y, nan)

# min_x = x.min()
# max_x = x.max()

# x = x.reshape(-1, 1)

# model.fit(x, y)
# # r_sq = model.score(x, y)
# # print(f"R: {r_sq}")
# # print(f"Slop: {model.coef_}")

# linspace = np.linspace(min_x, max_x)
# y_line = linspace * model.coef_ + model.intercept_ 

# import plotly.graph_objects as go

# fig = go.Figure()
# fig.add_trace(go.Scatter(x=x.flatten(), y=y, mode="markers"))
# # fig.add_trace(go.Scatter(x=linspace, y=y_line, mode="lines"))
# fig.show()






























# ################################################################################


# # import h5py
# # import os
# # import numpy as np
# # import pandas as pd
# # import numpyro as no
# # import plotly.express as px
# # import plotly.graph_objects as go

# # class Posterior():
# #     def __init__(self, chain_files, key, burnin, prob):
# #         ## chain_files: list of hdf5 objects
# #         ## key: key value for parameter
# #         ## burning: number of samples to exclude for burnin
# #         ## prob: hpd interval 
# #         self.param = key
# #         self.n_chains = len(chain_files)
# #         self.lower = round((1.0 - prob) / 2, 3) * 100 
# #         self.upper = 100 - self.lower
# #         if len(chain_files[0][key].shape) == 1:
# #             self.param_dims = 1
# #         elif len(chain_files[0][key].shape) == 2 or len(chain_files[0][key].shape) == 3:
# #             self.param_dims = 2 
# #         else:
# #             quit("Invalid input array shape")

# #         chains = []
# #         for chain_file in chain_files:
# #             if self.param_dims == 1:
# #                 chains.append(np.array(chain_file[key]))
# #             elif self.param_dims == 2:
# #                 chains.append(np.array(chain_file[key]).reshape(chain_file[key].shape[0], -1))
# #         self.chains = np.array(chains)

# #         data = self.chains[..., burnin:]
# #         if self.param_dims == 1:
# #             summary = no.diagnostics.summary(data, prob=prob)["Param:0"]
# #             summary.update({"parameter": "{}".format(key)})
# #             summary.move_to_end("parameter", last=False)
# #             self.summary = pd.DataFrame(summary, index=[0])
# #         elif self.param_dims == 2:
# #             summaries = []
# #             for i in range(data.shape[1]):
# #                 summary = no.diagnostics.summary(data[:,i,:], prob=prob)["Param:0"]
# #                 summary.update({"parameter": "{}_{}".format(key, i)})
# #                 summary.move_to_end("parameter", last=False)
# #                 summaries.append(summary)
# #             self.summary = pd.DataFrame(summaries)

# #     def to_df(self):
# #         data = self.chains
# #         df = pd.DataFrame(data.transpose())
# #         df.index.names = ["Sample"]
# #         df.columns = list(range(1, self.n_chains + 1)) 
# #         return df

# #     def to_long_df(self, indices=None):
# #         data = self.chains 
# #         if indices is not None:
# #             data = data[:,indices,:] 
# #         else:
# #             indices = range(data.shape[1])
# #         dfs = [] 
# #         for i in range(data.shape[1]):
# #             df = pd.DataFrame(data[:,i,:].transpose())
# #             df.columns = list(range(1, self.n_chains + 1)) 
# #             df.index.names = ["Sample"]
# #             df["param"] = indices[i]
# #             dfs.append(df)
# #         df = pd.concat(dfs)   
# #         return df

# #     def histogram(self, path, indices=None):
# #         if self.param_dims == 1:
# #             if indices is not None:
# #                 print("Warning: indices arg not used for {} parameter".format(self.param))
# #             df = self.to_df()
# #             fig = px.histogram(df, x=[*range(1, self.n_chains + 1)], title=self.param, 
# #                                barmode="overlay", histnorm="probability density")
# #         elif self.param_dims == 2:
# #             df = self.to_long_df(indices)
# #             fig = px.histogram(df, x=[*range(1, self.n_chains + 1)], title=self.param,
# #                     barmode="overlay", #histnorm='probability density',
# #                     facet_col="param", facet_col_wrap=4, facet_row_spacing=0.01, 
# #                     facet_col_spacing=0.01)
                     
# #         fig.update_layout(legend_title_text="Chain", title=dict(x=0.5))
# #         # fig.write_html(path)
# #         fig.show() 

# #     def trace(self, path, indices=None):
# #         if self.param_dims == 1:
# #             if indices is not None:
# #                 print("Warning: indices arg not used for {} parameter".format(self.param))
# #             df = self.to_df()
# #             fig = px.line(df, x=df.index, y=[*range(1, self.n_chains + 1)], 
# #                           title=self.param)
# #         elif self.param_dims == 2:
# #             df = self.to_long_df(indices)
# #             # figs = []
# #             # for name, group in df.groupby(by="param"):
# #             #     fig = px.line(
# #             #             group, 
# #             #             x=group.index, 
# #             #             y=[*range(1, self.n_chains + 1)], 
# #             #             height=200,
# #             #             width=200,
# #             #             title="{} {}".format(param, name))
# #             #     figs.append(fig)
# #             # with open(path, "a") as fh:
# #             #     for i in figs:
# #             #         fh.write(i.to_html(full_html=False, include_plotlyjs='cdn'))
                
# #             fig = px.line(df, x=df.index, y=[*range(1, self.n_chains + 1)], 
# #                     title=self.param, facet_col="param", facet_col_wrap=4, 
# #                     facet_row_spacing=0.01, facet_col_spacing=0.01)
# #         fig.update_layout(legend_title_text="Chain", title=dict(x=0.5))
# #         fig.update_traces(line=dict(width=1))
# #         # fig.write_html(path)
# #         fig.show()

# #     def covered(self, x):
# #         ## Returns list of booleans
# #         ## True if x is covered by credible interval
# #         ## False if x is not covered by credible interval
# #         lower = "{}%".format(self.lower)
# #         upper = "{}%".format(self.upper)
# #         df = self.summary
# #         return (x > df[lower]) & (x < df[upper])

# # def cline_func(a, b, h):
# #     thetas = h + 2 * h *(1 - h) * (a + b * (2 * h - 1))
# #     thetas[thetas > 1] = 1
# #     thetas[thetas < 0] = 0
# #     return thetas

# # def cline_plot(alpha, beta, alpha_outliers=None, beta_outliers=None):
# #     assert len(alpha) == len(beta)
# #     if alpha_outliers is None:
# #         alpha_outliers = pd.Series(False, index=range(len(alpha)))
# #     else:
# #         assert len(alpha) == len(alpha_outliers)
# #     if beta_outliers is None:
# #         beta_outliers = pd.Series(False, index=range(len(beta)))
# #     else:
# #         assert len(alpha) == len(beta_outliers)
# #     x = np.linspace(0, 1, 100) 
# #     fig = go.Figure()
# #     def add_trace( a, b, color):
# #         for i in range(len(a)): 
# #             phi = cline_func(a[i], b[i], x)  
# #             fig.add_trace(go.Scatter(x=x, y=phi, mode="lines", line=dict(color=color)))
# #     # Non outlier
# #     add_trace( 
# #         alpha[~alpha_outliers & ~beta_outliers].reset_index(drop=True),
# #         beta[~alpha_outliers & ~beta_outliers].reset_index(drop=True),
# #         "lightslategray")
# #     # alpha outlier
# #     add_trace(
# #         alpha[alpha_outliers & ~beta_outliers].reset_index(drop=True),
# #         beta[alpha_outliers & ~beta_outliers].reset_index(drop=True),
# #         "orange")
# #     # beta outlier
# #     add_trace(
# #         alpha[~alpha_outliers & beta_outliers].reset_index(drop=True),
# #         beta[~alpha_outliers & beta_outliers].reset_index(drop=True),
# #         "red")
# #     # alpha and beta outlier
# #     add_trace(
# #         alpha[alpha_outliers & beta_outliers].reset_index(drop=True),
# #         beta[alpha_outliers & beta_outliers].reset_index(drop=True),
# #         "green")
# #     fig.update_layout(
# #         showlegend=False,
# #         yaxis_range=[0,1],
# #         xaxis_range=[0,1],
# #         xaxis_title="Hybrid Index",
# #         yaxis_title="Prob. Ancestry",
# #         width=1000,
# #         height=1000,
# #     )
# #     fig.show()

# # from bgc import BGC, Posterior, cline_func, cline_plot





# # class Posterior():
# #     def __init__(self, paths, burnin=0, prob=0.9):
# #         chain_files = [h5py.File(path) for path in paths]
# #         self.keys = list(chain_files[0].keys())
# #         self.posterior_chains = {}
# #         self.n_chains = len(chain_files)
# #         self.burnin = burnin
# #         for key in self.keys:
# #             chains = []
# #             for chain_file in chain_files:
# #                 if len(chain_file[key].shape) == 1 or len(chain_file[key].shape) == 2:
# #                     chains.append(np.array(chain_file[key]))
# #                 elif len(chain_file[key].shape) == 3:
# #                     chains.append(np.array(chain_file[key]).reshape(chain_file[key].shape[0], -1))
# #                 else:
# #                     quit("Invalid input array shape")
# #             self.posterior_chains[key] = np.array(chains)
# #         self.get_summary(prob=prob)

#     # def get_summary(self, prob=0.9): 
#     #     self.prob = prob
#     #     self.summary = {}
#     #     for key in self.keys:
#     #         summaries = []
#     #         for key in self.keys:
#     #             data = self.posterior_chains[key][..., self.burnin:]
#     #             if len(data.shape) == 3:
#     #                 for i in range(data.shape[1]):
#     #                     summary = no.diagnostics.summary(data[:,i,:], prob=prob)["Param:0"]
#     #                     summary.update({"parameter": "{}_{}".format(key, i)})
#     #                     summary.move_to_end("parameter", last=False)
#     #                     summaries.append(summary)
#     #             else:
#     #                 summary = no.diagnostics.summary(data, prob=prob)["Param:0"]
#     #                 summary.update({"parameter": "{}".format(key)})
#     #                 summary.move_to_end("parameter", last=False)
#     #                 summaries.append(summary)
#     #         self.summary[key] = summaries



#     # def trace_plot(self, param, path, indices=None):
#     #     shape = self.posterior_chains[param].shape
#     #     if len(shape) == 3:
#     #         df = self.to_long_df(param, indices)
#     #         # figs = []
#     #         # for name, group in df.groupby(by="param"):
#     #         #     fig = px.line(
#     #         #             group, 
#     #         #             x=group.index, 
#     #         #             y=[*range(1, self.n_chains + 1)], 
#     #         #             height=200,
#     #         #             width=200,
#     #         #             title="{} {}".format(param, name))
#     #         #     figs.append(fig)
#     #         # with open(path, "a") as fh:
#     #         #     for i in figs:
#     #         #         fh.write(i.to_html(full_html=False, include_plotlyjs='cdn'))
                
#     #         fig = px.line(df, x=df.index, y=[*range(1, self.n_chains + 1)], 
#     #                 title=param, facet_col="param", facet_col_wrap=4, 
#     #                 facet_row_spacing=0.01, facet_col_spacing=0.01)
#     #     else:
#     #         if indices is not None:
#     #             print("Warning: indices arg not used for {} parameter".format(param))
#     #         df = self.to_df(param)
#     #         fig = px.line(df, x=df.index, y=[*range(1, self.n_chains + 1)], 
#     #                       title=param)
#     #     fig.update_layout(legend_title_text="Chain", title=dict(x=0.5))
#     #     fig.update_traces(line=dict(width=1))
#     #     fig.write_html(path)
#     #     # fig.show()

#     # def histogram(self, param, path, indices=None):
#     #     shape = self.posterior_chains[param].shape
#     #     if len(shape) == 3:
#     #         df = self.to_long_df(param, indices)
#     #         fig = px.histogram(df, x=[*range(1, self.n_chains + 1)], title=param,
#     #                 barmode="overlay", #histnorm='probability density',
#     #                 facet_col="param", facet_col_wrap=4, facet_row_spacing=0.01, 
#     #                 facet_col_spacing=0.01)
                     
#     #     else:
#     #         if indices is not None:
#     #             print("Warning: indices arg not used for {} parameter".format(param))
#     #         df = self.to_df(param)
#     #         fig = px.histogram(df, x=[*range(1, self.n_chains + 1)], title=param, 
#     #                            barmode="overlay", histnorm="probability density")
#     #     fig.update_layout(legend_title_text="Chain", title=dict(x=0.5))
#     #     fig.write_html(path)
#     #     # fig.show() 

#     # def summary_to_file(self, path):
#     #     # summaries = []
#     #     # for key in self.keys:
#     #     #     data = self.posterior_chains[key][..., self.burnin:]
#     #     #     if len(data.shape) == 3:
#     #     #         for i in range(data.shape[1]):
#     #     #             summary = no.diagnostics.summary(data[:,i,:], prob=prob)["Param:0"]
#     #     #             summary.update({"parameter": "{}_{}".format(key, i)})
#     #     #             summary.move_to_end("parameter", last=False)
#     #     #             summaries.append(summary)
#     #     #     else:
#     #     #         summary = no.diagnostics.summary(data, prob=prob)["Param:0"]
#     #     #         summary.update({"parameter": "{}".format(key)})
#     #     #         summary.move_to_end("parameter", last=False)
#     #     #         summaries.append(summary)
#     #     df = pd.DataFrame(self.summary) 
#     #     print(df)
#     #     # df.to_csv(path, index=False)
#     #     # print("Summary written to {}".format(path))


# # class Posterior():
# #     def __init__(self, name, chains):
# #         self.name = name
# #         self.chains = chains

# #     def get_param_array(self, param):
# #         chains = []
# #         for i in range(1, self.chains + 1): 
# #             fname = "{}_{}.hdf5".format(self.name, i) 
# #             f = h5py.File(os.path.join(self.name, fname), 'r')
# #             chains.append(np.array(f[param]).reshape(f[param].shape[0], -1))
# #         return np.array(chains)
 
# #     def to_df(self, data):
# #         df = pd.DataFrame(data.transpose())
# #         df.index.names = ["Sample"]
# #         return df

# #     def to_long_df(self, data):
# #         dfs = [] 
# #         for i in range(data.shape[1]):
# #             df = self.to_df(data[:,i,:]) 
# #             df["locus"] = i
# #             dfs.append(df)
# #         return pd.concat(dfs)
  
# #     def trace_plots(self, df, param): 
# #         fig = px.line(df, x=df.index, y=[*range(self.chains)],
# #                       title=param)
# #         fig.update_traces(line=dict(width=1))
# #         fig.update_layout(legend_title_text="Chain")
# #         fig.show()
    
# #     def trace_plots_facet(self, df, param):
# #         fig = px.line(df, x=df.index, y=[*range(self.chains)], facet_col="locus",
# #                       title=param, facet_col_wrap=4)
# #         fig.update_traces(line=dict(width=1))
# #         fig.update_layout(legend_title_text="Chain")
# #         fig.show()

# #     def histograms(self, df, param):
# #         fig = px.histogram(df, x=[*range(self.chains)], facet_col="locus", 
# #                            title=param, facet_col_wrap=4, 
# #                            histnorm="probability density", barmode="overlay")
# #         fig.update_layout(legend_title_text="Chain")
# #         fig.show()

# #     # def ess(self, data):
# #     #     values = []
# #     #     for i in range(data.shape[1]):
# #     #         values.append(no.diagnostics.effective_sample_size(data[:,i,:]))
# #     #     return values 

# #     # def rhat(self, data):
# #     #     values = []
# #     #     for i in range(data.shape[1]):
# #     #         values.append(no.diagnostics.gelman_rubin(data[:,i,:]))
    
# #     # def hpdi(self, data, prob):
# #     #     lower = []
# #     #     upper = []
# #     #     for i in range(data.shape[1]):
# #     #         hpdi = no.diagnostics.hpdi(data[:,i,:], prob)
# #     #         print(hpdi.shape)

# #     def summary(self, data, prob):
# #         summaries = []
# #         for i in range(data.shape[1]):
# #             summaries.append(no.diagnostics.summary(data[:,i,:])["Param:0"])
# #         df = pd.DataFrame(summaries)
# #         return df
    



# # ['LnL', 'alpha', 'beta', 'eta', 'eta-quantile', 'gamma-quantile', 
# # 'gamma-quantile-local', 'hi', 'interspecific-het', 'kappa', 'kappa-quantile', 
# # 'rho', 'tau-alpha', 'tau-beta', 'zeta-quantile', 'zeta-quantile-local']

# # f = h5py.File("clust-90-missing-0.9-exclude-1/clust-90-missing-0.9-exclude-1_1.hdf5", "r")
# # print(list(f.keys()))

# # P = Posterior("clust-90-missing-0.9-exclude-1", 5)

# # def get_param_array(param, n_chains, name):
# #     chains = []
# #     for i in range(1, n_chains + 1): 
# #         fname = "{}_{}.hdf5".format(name, i) 
# #         f = h5py.File(os.path.join(name, fname), 'r')
# #         chains.append(np.array(f[param]).reshape(f[param].shape[0], -1))
# #     return np.array(chains)







# # chain_files = [h5py.File(path) for path in paths]
# # param_keys = list(chain_files[0].keys())
# # posterior_chains = {} 

# # for key in param_keys:
# #     chains = []
# #     for chain_file in chain_files:
# #         if len(chain_file[key].shape) == 1 or len(chain_file[key].shape) == 2:
# #             chains.append(np.array(chain_file[key]))
# #         elif len(chain_file[key].shape) == 3:
# #             chains.append(np.array(chain_file[key]).reshape(chain_file[key].shape[0], -1))
# #         else:
# #             quit("Invalid input array shape")
# #     posterior_chains[key] = np.array(chains)








    
# # alpha = P.get_param_array("alpha")
# # beta = P.get_param_array("beta")
# # lnl = P.get_param_array("LnL")
# # eta = P.get_param_array( "eta")
# # eta_quantile = P.get_param_array( "eta-quantile")
# # gamma_quantile = P.get_param_array( "gamma-quantile")
# # gamma_quantile_local = P.get_param_array( "gamma-quantile-local")
# # hi = P.get_param_array( "hi")
# # interspecific_het = P.get_param_array( "interspecific-het")
# # kappa = P.get_param_array( "kappa")
# # kappa_quantile = P.get_param_array( "kappa-quantile")
# # rho = P.get_param_array( "rho")
# # tau_alpha = P.get_param_array( "tau-alpha")
# # tau_beta = P.get_param_array( "tau-beta")
# # zeta_quantile = P.get_param_array( "zeta-quantile")
# # zeta_quantile_local = P.get_param_array( "zeta-quantile-local")

# # print(P.to_df(lnl[0,:,:]))

# # P.histograms(P.to_long_df(lnl), "LnL")



# # P.histograms(P.to_long_df(alpha[:,0:20,:]), "Alpha")
# # P.histograms(P.to_long_df(beta[:,0:20,:]), "Beta")
# # P.trace_plots_facet(P.to_long_df(alpha[:,0:20,:]), "Alpha")
# # P.trace_plots_facet(P.to_long_df(beta[:,0:20,:]), "Beta")

# # alpha_df = P.summary(alpha, 0.95)
# # beta_df = P.summary(beta, 0.95)
# # alpha_df.to_csv("summary-alpha.csv")
# # beta_df.to_csv("summary-beta.csv")





# # print(beta_df.columns)

# # x = np.linspace(0, 1, 100) 
# # # theta = cline_func(1, -1, x)

# # fig = go.Figure()
# # for i in range(len(alpha_df)): 
# #     phi = cline_func(alpha_df.iloc[i]["mean"], beta_df.iloc[i]["mean"], x)  
# #     fig.add_trace(go.Scatter(x=x, y=phi, mode="lines", line=dict(color="black")))

# # fig.update_layout(
# #     yaxis_range=[0,1],
# #     xaxis_range=[0,1])
# # # fig.layout.yaxis.scaleanchor = "x"
# # fig.show()

