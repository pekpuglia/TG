using Printf
function transfer_table(transfer::Transfer, tspan_ppdot::Tuple, L, T, tolepsilon, scenario)
    format = raw"
\begin{table}[htpb]
    \centering
    \begin{tabular}{cccc} \toprule
    \multicolumn{2}{c}{\textbf{Maneuver type}} & \multicolumn{2}{c}{\texttt{MANEUVER_TYPE}} \\ \midrule
    \(L\) (m) & \(T\) (s) & \(\varepsilon\) & \(\lVert \Delta \pos_{f} \rVert\) (m)    \\ \midrule
    SCALEL          & SCALET          & TOLEPSILON                & DELTAXFINAL                        \\ \midrule
    \(\max \lVert p \rVert\) & MAXNORMP     & \textbf{Diagnostic}   & DIAGNOSTIC        \\ \midrule
    \textbf{Impulse} & \(t\) (s) & \(\Delta v\) (m/s) & \(1 - p \cdot \hat{u}\) \\ \midrule
    IMPINDEX                 & IMPTIME          & IMPDV             & IMPDIR                    \\
    \textbf{Total}   & TOTALT          & TOTALDV             &                     \\ \bottomrule   
    \end{tabular}
    \caption{Summary of optimization for \texttt{MANEUVER_TYPE} MODEL SCENARIO.}
    \label{tab:tb_SCEN_MANEUVER_TYPE_tab}
\end{table}
"

    deltaxfinal = let 
        recomp_transfer = recompute_model(transfer, RK8)
        norm(coasts(recomp_transfer)[end].rcoast[:, end] - transfer.X2[1:3])
    end
    p = tspan_ppdot[2][1:3, :]
    maxnormp = maximum(norm.(eachcol(p)))
    diagnostic = diagnose_ppdot(tspan_ppdot)
    imp_ts = impulse_times(transfer)
    total_t = transfer.transfer_time
    totaldv = total_dV(transfer)
    
    partial_table = replace(format, 
        "MANEUVER_TYPE" => transfer_type(transfer),
        "SCALEL"      => round(L, digits=3),
        "SCALET"      => round(T, digits=3),
        "TOLEPSILON"  => @sprintf("%.2e", tolepsilon),
        "DELTAXFINAL" => @sprintf("%.5e", deltaxfinal),
        "MAXNORMP"    => round(maxnormp, sigdigits=5),
        "DIAGNOSTIC"  => diagnostic,
        "TOTALT"      => round(total_t, digits=5),
        "TOTALDV"     => round(totaldv, digits=5),
        "SCENARIO"    => scenario,
        "SCEN"        => (scenario |> lowercase |> split) .|> first |> join,
        "MODEL"       => model_name(transfer.model)
    )
    
    tspan = tspan_ppdot[1]
    #assumes impulse times are always on tspan - check
    index_p_on_imp = [findfirst(==(imp_t), tspan) for imp_t = imp_ts]
    p_on_imp = p[:, index_p_on_imp] |> eachcol |> collect
    #add impulses to table
    partial_table_lines = split(partial_table, "\n")
    impulse_lines = repeat([partial_table_lines[10]], transfer.nimp)
    for (i, t, imp, p) = zip(1:transfer.nimp, imp_ts, impulses(transfer), p_on_imp)
        dv = imp.deltaVmag
        impulse_lines[i] = replace(impulse_lines[i], 
            "IMPINDEX" => i, 
            "IMPTIME" => round(t, digits=5), 
            "IMPDV" => round(dv, digits=5), 
            "IMPDIR" => round(1 - dot(p, imp.deltaVdir), digits=3)
        )
    end
    impulse_lines[end] = impulse_lines[end]*"\\midrule"
    final_table = join([partial_table_lines[1:9]..., impulse_lines..., partial_table_lines[11:end]...], "\n")
    final_table
end

function setup_table(start_orb::KeplerianElements, final_orb::KeplerianElements, transfer_time)
    format = raw"
\begin{table}[htbp]
    \centering
    \begin{tabular}{ccc} \toprule
        Element & Initial & Final \\ \midrule
        \(a\)      & \(A1\) km         & \(A2\) km   \\
        \(e\)      & \(E1\)            & \(E2\)        \\
        \(i\)      & \(I1^\circ\)      & \(I2^\circ\) \\
        \(\Omega\) & \(RAAN1^\circ\)   & \(RAAN2^\circ\)  \\
        \(\omega\) & \(OMEGA1^\circ\)  & \(OMEGA2^\circ\)  \\
        \(\theta\) & \(THETA1^\circ\)  & \(THETA2^\circ\)  \\ 
        Transfer time & \multicolumn{2}{c}{\(TRANSFER_TIME\)} \\\bottomrule
    \end{tabular}
    \caption{Orbital elements used for the TODO transfer case analysis}
    \label{tab:TODO_orb_elems}
\end{table}
"
    replace(format,
    "A1"     => round(start_orb.a / 1e3, digits = 3),
    "E1"     => round(start_orb.e, digits = 3),
    "I1"     => round(rad2deg(start_orb.i), digits = 3),
    "RAAN1"  => round(rad2deg(start_orb.Ω), digits = 3),
    "OMEGA1" => round(rad2deg(start_orb.ω), digits = 3),
    "THETA1" => round(rad2deg(start_orb.f), digits = 3),
    "A2"     => round(final_orb.a / 1e3, digits = 3),
    "E2"     => round(final_orb.e, digits = 3),
    "I2"     => round(rad2deg(final_orb.i), digits = 3),
    "RAAN2"  => round(rad2deg(final_orb.Ω), digits = 3),
    "OMEGA2" => round(rad2deg(final_orb.ω), digits = 3),
    "THETA2" => round(rad2deg(final_orb.f), digits = 3),
    "TRANSFER_TIME" => round(transfer_time, digits = 3)
    )
end

function summary_table(transfers::Transfer...)
    
end

function dump_table(filename::String, dlm::String, table::String)
    file_contents = read(filename, String)

    split_file = split(file_contents, dlm)
    length(split_file) |> display
    new_file = split_file[1]*dlm*table*dlm*split_file[3]

    open(filename, "w") do f
        write(f, new_file)
    end
end