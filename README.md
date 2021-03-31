# fibra
FIBRA - Fixed Income Brazil. Government and Corporate Bonds Pricing.

O módulo permite precificar títulos públicos federais e corporativos do Brasil.

Com o Fibra é possível calcular o preço dos seguintes títulos | Bonds available:

### Títulos do Governo | Government Bonds
- NTN-B: Nota do Tesouro Nacional Série B.
- NTN-F: Nota do Tesouro Nacional Série F.
- LTN: Letra do Tesouro Nancional.

### Títulos Corporativos - Debêntures | Corporate Bonds
- Percentual CDI.
- DI Spread.
- IPCA Spread.

No jupyter notebook comece por | Start with jupyter notebook:

`%run ../docs/fibra.py`

## Notas do Tesouro Nacional | Bonds

`bond = Bond(date="2021-03-29", maturity="2050-08-15", ytm=4.3023, coupon=6, freq=2)`

### Calendário de pagamentos | Settlement schedule

`bond.schedule()`

| Date       |
| -----------|
|	2021-08-16 |
|	2022-02-15 |
|	2022-08-15 |
|	2023-02-15 |
|	2023-08-15 |
|	2024-02-15 |
|	2024-08-15 |
| 2025-02-17 |
| 2025-08-15 |
| 2026-02-18 |

### Fluxo de caixa | Cash flow

`bond.cashflow()`

|       Date | Days |    Years | Cash Flows | PV of Cash Flows |
|-----------:|-----:|---------:|-----------:|-----------------:|
| 2021-08-16 |   97 | 0.384921 |     2.9563 |         2.908753 |
| 2022-02-15 |  224 | 0.888889 |     2.9563 |         2.847654 |
| 2022-08-15 |  348 | 1.380952 |     2.9563 |         2.789238 |
| 2023-02-15 |  476 | 1.888889 |     2.9563 |         2.730193 |
| 2023-08-15 |  599 | 2.376984 |     2.9563 |         2.674633 |
| 2024-02-15 |  723 | 2.869048 |     2.9563 |         2.619766 |
| 2024-08-15 |  850 | 3.373016 |     2.9563 |         2.564738 |
| 2025-02-17 |  979 | 3.884921 |     2.9563 |         2.510026 |
| 2025-08-15 | 1102 | 4.373016 |     2.9563 |         2.458946 |

### Dados do título | Bond overview

`bond.des()`

|                     |        Days |
|--------------------:|------------:|
|          Convention |     BUS/252 |
|           Frequency | Semi-annual |
|        Current date |  2021-03-29 |
|       Maturity date |  2050-08-15 |
|    Last coupon date |  2021-02-15 |
|    Next coupon date |  2021-08-16 |
| Days to next coupon |          97 |
|  Remaining payments |          59 |
|    Term to maturity |       29.29 |
|         Coupon rate |        6.00 |
|        Coupon value |      2.9563 |
|       Nominal value |      100.00 |
|               Price |  128.544273 |
|   Yield to matutity |      4.3023 |
|       Current yield |      4.6677 |
|            Duration |     15.8254 |
|   Modified Duration |     15.1726 |
|                DV01 |       0.195 |
|           Convexity |    341.0802 |
|         Sensibility |    -13.4672 |

## Debêntures | Corporate Bonds

### Parâmetros | Required parameters

- ***DATA***: Data de referência para o cálculo do preço do ativo.
- ***VNE***: Valor Nominal de Emissão.
- ***VNA***: Valor Nominal Atualizado na data de referência.
- ***PU***: Preço Unitário do título ao Par (marcado na curva).
- ***TAXA***: Taxa indicativa (taxa de mercado) da ANBIMA.
- ***FREQ***: Frequência de pagamentos (1: Anual, 2: Semestral, 3: Trimestral, etc).

DATA_natura = DATA_REF<br/>
VNE_natura = 10000
VNA_natura = 10000<br/>
PU_natura = 10170.808970
TAXA_natura = 0.7883<br/>
FREQ_natura = 2<br/>

Calendário de eventos | Events calendar:  
`cal_natura = {"2021-09-27": 1}`

`natura = DebentureSpread(date=DATA_natura, maturity="2021-09-25", vne=VNE_natura, vna=VNA_natura, pu=PU_natura,  
                          issue_spread=1.75, market_spread=TAXA_natura, freq=FREQ_natura, redemption=cal_natura,  
                          yield_curve_file=YIELD_CURVE_PATH)`
                 
### Dados do título | Bond overview

`natura.des()`
|                      | Debenture DI Spread |
|---------------------:|--------------------:|
|           Convention |           DI Spread |
|            Frequency |         Semi-annual |
|                 Date |          2021-03-18 |
|             Maturity |          2021-09-25 |
|    Last payment date |          2020-09-25 |
|    Next payment date |          2021-03-25 |
| Days to next payment |                   5 |
|   Remaining payments |                   2 |
|     Term to maturity |                0.53 |
|                  VNE |        10000.000000 |
|                  VNA |        10000.000000 |
|               PU Par |        10170.808970 |
|                Price |        10221.045950 |
|      Price/Par ratio |             100.49% |
|           Issue rate |             1.7500% |
|          Market rate |             0.7883% |
|             Duration |                0.52 |

### Fluxo de caixa | Cash flow

`natura.cashflow()`
|       Date | Days |       DI |     VNA |  Interests | Redemptions |     Payments | PV of Cash Flows |
|-----------:|-----:|---------:|--------:|-----------:|------------:|-------------:|-----------------:|
| 2021-03-25 |    5 | 2.645583 | 10000.0 | 179.583171 |         0.0 |   179.583171 |       179.462193 |
| 2021-09-27 |  133 | 3.974169 | 10000.0 | 292.831916 |     10000.0 | 10292.831916 |     10041.583758 |

Veja alguns exemplos de uso em Bonds.ipynb e Debentures.ipynb.

See more in Bonds.ipynb and Debentures.ipynb.
