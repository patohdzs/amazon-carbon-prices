import numpy as np 
import pandas as pd
import os
from hmmlearn import hmm
from scipy.linalg import logm, expm
import warnings
import contextlib
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

@contextlib.contextmanager
def suppress_output():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr




def format_latex_table(result_uncon, result_con):
    _, _, _, mus_uncon, sigmas_uncon, _, _, _, P_uncon = result_uncon
    _, _, _, mus_con, sigmas_con, _, _, _, P_con = result_con

    table = r"""\begin{tabular}{cccccc}
\toprule
&\multicolumn{2}{c}{distinct variances} & &\multicolumn{2}{c}{a common variance} \\
\midrule
& low price & high price & & low price & high price \\
"""
    table += f" & {mus_uncon[0]:.2f} & {mus_uncon[1]:.2f} & & {mus_con[0]:.2f} & {mus_con[1]:.2f} \\\\\n"
    table += f"s.d. & {sigmas_uncon[0]:.3f} & {sigmas_uncon[1]:.3f} & & {sigmas_con[0]:.3f} & {sigmas_con[1]:.3f} \\\\\n"
    table += r"""\midrule
& & & transition probabilities & & \\
& low & high & & low & high \\
"""
    table += f"low  & {P_uncon[0,0]:.3f} & {P_uncon[0,1]:.3f} & & {P_con[0,0]:.3f} & {P_con[0,1]:.3f} \\\\\n"
    table += f"high & {P_uncon[1,0]:.3f} & {P_uncon[1,1]:.3f} & & {P_con[1,0]:.3f} & {P_con[1,1]:.3f} \\\\\n"
    table += r"\bottomrule" + "\n\\end{tabular}"

    return table

def information_criteria(result_uncon, result_con):
    aic_u, ll_u, bic_u = result_uncon[:3]
    aic_c, ll_c, bic_c = result_con[:3]

    table = r"""\begin{tabular}{ccc}
\hline
 & distinct variances & common variance \\
\hline
"""
    table += f"log likelihood & {ll_u:.2f} & {ll_c:.2f} \\\\\n"
    table += f"aic & {aic_u:.2f} & {aic_c:.2f} \\\\\n"
    table += f"bic & {bic_u:.2f} & {bic_c:.2f} \\\\\n"
    table += r"\hline" + "\n\\end{tabular}"

    return table


def plot_hmm_results(mu,predict, price, var):
    
    if mu[0]>mu[1]:
        predict = predict[:,0]
    else:
        predict = predict[:,1]
    predict = predict 
    price = price  

    start_date = '1995-01'
    dates = pd.date_range(start=start_date, periods=276, freq='M')

    df = pd.DataFrame({
    'date': dates,
    'predict': predict
    })

    # Save to CSV
    df.to_csv(f'output/tables/smooth_prob_{var}.csv', index=False)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the first dataset with the left y-axis (ax1)
    ax1.plot(dates, predict, linewidth=4, label='high state probability', color='b')
    ax1.set_xlabel('time')
    ax1.set_ylabel('smoothed probability', color='b',fontsize=16)
    ax1.tick_params(axis='y', labelcolor='b')

    ax1.xaxis.set_major_locator(mdates.YearLocator(base=5))  # Major ticks every 5 years
    ax1.xaxis.set_minor_locator(mdates.YearLocator(base=1))  # Minor ticks every year
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    # Create a second y-axis on the right
    ax2 = ax1.twinx()

    # Plot the second dataset with the right y-axis (ax2)
    ax2.plot(dates, price, linewidth=4, label='Actual Price', color='r')
    ax2.set_ylabel('actual Price', color='r',fontsize=16)
    ax2.tick_params(axis='y', labelcolor='r')

    if mu[0]>mu[1]:
        ax2.axhline(y=mu[0], color='g', linestyle='--', label='high price')
        ax2.axhline(y=mu[1], color='orange', linestyle='--', label='low price')
    else:
        ax2.axhline(y=mu[0], color='orange', linestyle='--', label='low price')
        ax2.axhline(y=mu[1], color='g', linestyle='--', label='high price')

    # Add legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines = lines1 + lines2
    labels = labels1 + labels2
    ax1.legend(lines, labels, loc='lower left',fontsize=13)

    plt.xticks(rotation=45)
    plt.savefig(f'output/figures/smooth_prob_{var}.png', format='png', bbox_inches='tight')
    plt.close()

def est(s_low=0.5,s_high=0.5,var='uncon'):
    
    warnings.filterwarnings("ignore", message="KMeans is known to have a memory leak on Windows with MKL, when there are less chunks than available threads. You can avoid it by setting the environment variable OMP_NUM_THREADS=2.")

    
    # need to input initial state prob s_low and s_high
    data_folder = os.getcwd()+"/data/calibration/"
    
    # read data
    df = pd.read_csv(data_folder+"seriesPriceCattle_prepared.csv")
    price = df['price_real_mon_cattle'].values.astype(float) 
    logprice=np.log(price)
    
    
    # estimating the model
    
    Q=logprice.reshape(logprice.shape[0],1)

    np.random.seed(12345)

    best_score = best_model = None
    n_fits = 500

    for idx in range(n_fits):

        if var=='uncon':
            model = hmm.GaussianHMM(n_components=2,random_state=idx,init_params='tmc',params='stmc')
        else:
            model = hmm.GaussianHMM(n_components=2,random_state=idx,init_params='tmc',params='stmc',covariance_type='tied')
            
        model.startprob_= np.array([s_low, s_high])
        with suppress_output():
            model.fit(Q)

        score = model.score(Q)

        if best_score is None or score > best_score:
            best_model = model
            best_score = score
      
    aic = best_model.aic(Q)
    bic = best_model.bic(Q)
    mus = np.exp(np.ravel(best_model.means_))  
    sigmas = np.ravel(np.sqrt([np.diag(c) for c in best_model.covars_]))
    P = np.round(best_model.transmat_,3)
    
    sorted_indices = np.argsort(mus)
    mus_sorted = mus[sorted_indices]
    sigmas_sorted = sigmas[sorted_indices]
    P_sorted = P[sorted_indices][:, sorted_indices]
    predict=best_model.predict_proba(Q)
    if var=='uncon':
        ll = (best_model.aic(Q)-12)/(-2)
    else:
        ll = (best_model.aic(Q)-10)/(-2)
        
        
    plot_hmm_results(mus,predict, price, var)
        
    
    return(aic,ll,bic,mus_sorted,sigmas_sorted,P_sorted)




def stationary(transition_matrix,mus):
    
    eigenvals, eigenvects = np.linalg.eig(transition_matrix.T)
    close_to_1_idx = np.isclose(eigenvals,1)
    target_eigenvect = eigenvects[:,close_to_1_idx]
    target_eigenvect = target_eigenvect[:,0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect)
    stationary_price=mus[0]*stationary_distrib[0]+mus[1]*stationary_distrib[1] 
    return(stationary_distrib,stationary_price)


def annual_transition(transition_matrix):
    dτ = 1/12  # time step
    P = transition_matrix  #probability transition matrix
    M = logm(P)/dτ  # instantenous generator
    P = expm(M)  # Updated Probability transition matrix
    return P


def iteration_est(initial_prob, num_iterations=5,var='uncon'):
    
    
    for i in range(num_iterations):
        if i == 0:
            s_low, s_high = initial_prob[0], initial_prob[1]
        else:
            s_low, s_high = sta_dist[0], sta_dist[1]
        
        aic, ll, bic, mus, sigmas, P = est(s_low, s_high, var)
        sta_dist, sta_price = stationary(P, mus)
        annual_P=annual_transition(P)
        
    return aic,ll,bic,mus,sigmas,P,sta_dist,sta_price,annual_P
    
    



result_uncon = iteration_est([0.5, 0.5], num_iterations=5, var='uncon')
result_con = iteration_est([0.5, 0.5], num_iterations=5, var='con')

print("distinct variances",   result_uncon)
print("common variances",   result_con)

latex_table = format_latex_table(result_uncon, result_con)
information_table = information_criteria(result_uncon, result_con)


# Save to file
with open("output/tables/hmm_results_table.tex", "w") as f:
    f.write(latex_table)
    

with open("output/tables/hmm_information_criteria.tex", "w") as f:
    f.write(information_table)
