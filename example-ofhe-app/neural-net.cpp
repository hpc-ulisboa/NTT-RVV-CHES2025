#define PROFILE

#include "openfhe.h"
#include "../ntt/ntt-rvv.h"

using namespace lbcrypto;


int neuralnetwork(uint32_t numSlots)
{
    // uint32_t numSlots = 1024;
    CCParams<CryptoContextCKKSRNS> parameters;

    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(1 << 10);

    parameters.SetNumLargeDigits(3);
    parameters.SetKeySwitchTechnique(HYBRID);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;
    usint dcrtBits = 59;
    usint firstMod = 60;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    std::vector<uint32_t> levelBudget = {1, 1};
    std::vector<uint32_t> bsgsDim = {0, 0};

    uint32_t levelsAvailableAfterBootstrap = 10;
    usint depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    // Generate crypto context.
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    // Enable features that you wish to use. Note, we must enable FHE to use bootstrapping.
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    usint ringDim = cryptoContext->GetRingDimension();
    std::cout << "CKKS scheme is using ring dimension " << ringDim << std::endl
              << std::endl;

    cryptoContext->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots);

    printf("bootstrap evaluated\n");

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);

    printf("keys generated\n");

    uint32_t first_layer_dim = 2;
    uint32_t second_layer_dim = 2;
    uint32_t third_layer_dim = 2;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.05, 0.05);

    std::vector<double> x;
    for (size_t i = 0; i < numSlots; i++)
    {
        x.push_back(dis(gen));
    }

    printf("input generated\n");

    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x, 1, depth - 1, nullptr, numSlots);
    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
    auto ciphertextAfter = cryptoContext->EvalBootstrap(ciph);

    printf("input encrypted\n");

    std::vector<std::vector<double>> weights_firstlayer;
    for (uint32_t j = 0; j < first_layer_dim; j++)
    {
        std::vector<double> weights1;
        for (size_t i = 0; i < numSlots; i++)
        {
            weights1.push_back(dis(gen));
        }
        weights_firstlayer.push_back(weights1);
    }

    std::vector<std::vector<double>> weights_secondlayer;
    for (uint32_t j = 0; j < second_layer_dim; j++)
    {
        std::vector<double> weights2;
        for (size_t i = 0; i < first_layer_dim; i++)
        {
            weights2.push_back(dis(gen));
        }
        weights_secondlayer.push_back(weights2);
    }

    std::vector<std::vector<double>> weights_thirdlayer;
    for (uint32_t j = 0; j < third_layer_dim; j++)
    {
        std::vector<double> weights3;
        for (size_t i = 0; i < second_layer_dim; i++)
        {
            weights3.push_back(dis(gen));
        }
        weights_thirdlayer.push_back(weights3);
    }

    printf("matrices generated\n");

    std::vector<Ciphertext<DCRTPoly>> ctweights;
    for (uint32_t i = 0; i < first_layer_dim; ++i)
    {

        Plaintext pweights = cryptoContext->MakeCKKSPackedPlaintext(weights_firstlayer[i], 1, depth - 1, nullptr, numSlots);

        ctweights.push_back(cryptoContext->Encrypt(keyPair.publicKey, pweights));
        ctweights[i] = cryptoContext->EvalBootstrap(ctweights[i]);
    }
    printf("weights encrypted, entering measured code\n");

    uint64_t start_cycles, end_cycles;
    __asm__ __volatile__("rdcycle %0" : "=r"(start_cycles));

    std::vector<Ciphertext<DCRTPoly>> ct_output_first_layer;
    for (uint32_t i = 0; i < first_layer_dim; i++)
    {
        auto finalResult = cryptoContext->EvalInnerProduct(ciphertextAfter, ctweights[i], numSlots);

        // Use Chebyshev as an activation function
        ct_output_first_layer.push_back(cryptoContext->EvalChebyshevFunction([](double x) -> double
                                                                             { if (x < 0) return 0; else return 1 * x; }, finalResult, -5, 5, 119));
        ct_output_first_layer[i] = cryptoContext->EvalBootstrap(ct_output_first_layer[i]);
    }

    std::vector<Ciphertext<DCRTPoly>> ct_output_second_layer;
    for (uint32_t i = 0; i < second_layer_dim; i++)
    {
        auto ct_output_second_layer_aux = cryptoContext->EvalMult(ct_output_first_layer[0], weights_secondlayer[i][0]);
        for (uint32_t j = 1; j < first_layer_dim; j++)
        {
            auto aux_ct = cryptoContext->EvalMult(ct_output_first_layer[j], weights_secondlayer[i][j]);
            ct_output_second_layer_aux = cryptoContext->EvalAdd(ct_output_second_layer_aux, aux_ct);
            // Use Chebyshev as an activation function
        }
        ct_output_second_layer_aux = cryptoContext->EvalChebyshevFunction([](double x) -> double
                                                                          { if (x < 0) return 0; else return 1 * x; }, ct_output_second_layer_aux, -5, 5, 119);
        ct_output_second_layer.push_back(cryptoContext->EvalBootstrap(ct_output_second_layer_aux));
    }

    std::vector<Ciphertext<DCRTPoly>> ct_output_third_layer;
    for (uint32_t i = 0; i < third_layer_dim; i++)
    {
        auto ct_output_third_layer_aux = cryptoContext->EvalMult(ct_output_second_layer[0], weights_thirdlayer[i][0]);
        for (uint32_t j = 1; j < second_layer_dim; j++)
        {
            auto aux_ct = cryptoContext->EvalMult(ct_output_second_layer[j], weights_thirdlayer[i][j]);
            ct_output_third_layer_aux = cryptoContext->EvalAdd(ct_output_third_layer_aux, aux_ct);
            // Use Chebyshev as an activation function
        }
        ct_output_third_layer_aux = cryptoContext->EvalChebyshevFunction(
            [](double x) -> double
            { if (x < 0) return 0; else return 1 * x; }, ct_output_third_layer_aux, -5, 5, 119);
        ct_output_third_layer.push_back(cryptoContext->EvalBootstrap(ct_output_third_layer_aux));
    }


    __asm__ __volatile__("rdcycle %0" : "=r"(end_cycles));
    std::cout << "Cycles: " << end_cycles - start_cycles << std::endl;

    // third layer results
    for (uint32_t i = 0; i < ct_output_third_layer.size(); ++i)
    {
        lbcrypto::Plaintext result;
        cryptoContext->Decrypt(keyPair.secretKey, ct_output_third_layer[i], &result);
        result->SetLength(numSlots);
    }

    return 0;
}


int main()
{
    neuralnetwork(1);

    free_ntts_mem();
    return 0;
}
