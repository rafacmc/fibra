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

### Fluxo de caixa | Cash flow

`bond.cashflow()`

### Dados do título | Bond overview

`bond.des()`

## Debêntures | Corporate Bonds

### Parâmetros | Required parameters

- ***DATA***: Data de referência para o cálculo do preço do ativo.
- ***VNE***: Valor Nominal de Emissão.
- ***VNA***: Valor Nominal Atualizado na data de referência.
- ***PU***: Preço Unitário do título ao Par (marcado na curva).
- ***TAXA***: Taxa indicativa (taxa de mercado) da ANBIMA.
- ***FREQ***: Frequência de pagamentos (1: Anual, 2: Semestral, 3: Trimestral, etc).

DATA_natura = DATA_REF

VNE_natura = 10000

VNA_natura = 10000

PU_natura = 10170.808970

TAXA_natura = 0.7883

FREQ_natura = 2

Calendário de eventos | Events calendar: `cal_natura = {"2021-09-27": 1}`

`natura = DebentureSpread(date=DATA_natura, maturity="2021-09-25", vne=VNE_natura, vna=VNA_natura, pu=PU_natura,
                          issue_spread=1.75, market_spread=TAXA_natura, freq=FREQ_natura, redemption=cal_natura,
                          yield_curve_file=YIELD_CURVE_PATH)`
                 
### Calendário de pagamentos | Settlement schedule

`natura.schedule()`

### Fluxo de caixa | Cash flow

`natura.cashflow()`

### Dados do título | Bond overview

`natura.des()`

Veja alguns exemplos de uso em `Bonds.ipynb` e `Debentures.ipynb`. See more in Bonds.ipynb` and `Debentures.ipynb`.
