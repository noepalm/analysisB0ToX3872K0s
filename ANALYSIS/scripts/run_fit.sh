for channel in Psi2S; do
    ./MVAcutter /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_${channel}.root -2 0.5 ${channel} 2022
done
# for channel in X3872; do
#     ./MVAcutter /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_${channel}.root -2 0.5 ${channel} 2022
# done

# ./MVAoptimization /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022_BDTtraining_CV3_X3872.root -4 0.4 0.1 0.05