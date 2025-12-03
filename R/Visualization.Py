"""
Módulo de visualización
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict

def plot_3_strategies_comparison(trading_results: Dict, train_start: int, 
                                test_start: int, test_end: int) -> None:
    """
    Genera gráficas para comparar 3 estrategias
    
    Parámetros:
    -----------
    trading_results : dict
        Resultados de trading
    train_start : int
        Inicio del train
    test_start : int
        Inicio del test
    test_end : int
        Fin del test
    """
    print("\nGENERANDO GRÁFICAS PARA LAS 3 ESTRATEGIAS...")
    print("-" * 50)
    
    # Gráfica 1: Spreads con movimientos y señales
    plt.figure(figsize=(15, 12))
    
    # Subgráfica 1: Spreads con señales
    plt.subplot(2, 1, 1)
    
    colores = {'beta_original': 'blue', 'beta_johansen': 'green', 'beta_largevar': 'red'}
    estilos = {'beta_original': '-', 'beta_johansen': '--', 'beta_largevar': '-.'}
    
    for metodo, resultado in trading_results.items():
        tiempo_test = np.arange(test_start, test_end)
        
        # Graficar spread
        plt.plot(tiempo_test, resultado['spread_test'],
                color=colores[metodo], linestyle=estilos[metodo],
                linewidth=2, label=f'{metodo}')
        
        # Marcar señales
        long_signals = np.where(resultado['position'] == 1)[0]
        short_signals = np.where(resultado['position'] == -1)[0]
        
        if len(long_signals) > 0:
            plt.scatter(tiempo_test[long_signals], resultado['spread_test'][long_signals],
                       color='green', s=100, marker='^',
                       label=f'{metodo} - Long' if metodo == list(trading_results.keys())[0] else "")
        
        if len(short_signals) > 0:
            plt.scatter(tiempo_test[short_signals], resultado['spread_test'][short_signals],
                       color='red', s=100, marker='v',
                       label=f'{metodo} - Short' if metodo == list(trading_results.keys())[0] else "")
    
    plt.axvline(x=test_start, color='black', linestyle='--', alpha=0.7, label='Inicio Trading')
    plt.title('SPREADS Y SEÑALES DE TRADING - 3 ESTRATEGIAS (T=200)', fontsize=14, fontweight='bold')
    plt.xlabel('Tiempo', fontsize=12)
    plt.ylabel('Valor del Spread', fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Subgráfica 2: Z-scores
    plt.subplot(2, 1, 2)
    
    for metodo, resultado in trading_results.items():
        tiempo_test = np.arange(test_start, test_end)
        plt.plot(tiempo_test, resultado['z_score'],
                color=colores[metodo], linestyle=estilos[metodo],
                linewidth=1.5, label=f'{metodo}', alpha=0.8)
    
    # Líneas de umbral
    plt.axhline(y=2, color='red', linestyle='--', alpha=0.7, label='Umbral Short (+2)')
    plt.axhline(y=1.5, color='orange', linestyle='--', alpha=0.5, label='Umbral Short (+1.5)')
    plt.axhline(y=-1.5, color='lightgreen', linestyle='--', alpha=0.5, label='Umbral Long (-1.5)')
    plt.axhline(y=-2, color='green', linestyle='--', alpha=0.7, label='Umbral Long (-2)')
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    plt.axvline(x=test_start, color='black', linestyle='--', alpha=0.7)
    
    plt.title('Z-SCORES Y UMBRALES DE TRADING - 3 ESTRATEGIAS', fontsize=14, fontweight='bold')
    plt.xlabel('Tiempo', fontsize=12)
    plt.ylabel('Z-Score', fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('spreads_3_estrategias_T200.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Gráfica 2: Retorno total
    plt.figure(figsize=(15, 10))
    
    # Subgráfica 1: Retornos acumulados
    plt.subplot(2, 2, 1)
    
    for metodo, resultado in trading_results.items():
        tiempo_retornos = np.arange(test_start, test_start + len(resultado['cumulative_returns']))
        plt.plot(tiempo_retornos, resultado['cumulative_returns'],
                color=colores[metodo], linewidth=3,
                label=f'{metodo} - Retorno: {resultado["total_return"]:.4f}')
    
    plt.title('RETORNOS ACUMULADOS - 3 ESTRATEGIAS', fontsize=14, fontweight='bold')
    plt.xlabel('Tiempo', fontsize=12)
    plt.ylabel('Retorno Acumulado', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Subgráfica 2: Comparación de métricas
    plt.subplot(2, 2, 2)
    
    metodos_estrategias = list(trading_results.keys())
    metricas = ['Retorno Total', 'Sharpe Ratio', 'Win Rate (%)']
    
    x_pos = np.arange(len(metricas))
    ancho_barra = 0.25
    
    for i, metodo in enumerate(metodos_estrategias):
        if metodo in trading_results:
            resultado = trading_results[metodo]
            valores = [
                resultado['total_return'],
                resultado['sharpe_ratio'],
                resultado['win_rate'] * 100
            ]
            plt.bar(x_pos + i * ancho_barra, valores, ancho_barra,
                   color=colores[metodo], alpha=0.7, label=metodo)
    
    plt.xticks(x_pos + ancho_barra, metricas)
    plt.title('COMPARACIÓN DE MÉTRICAS', fontsize=14, fontweight='bold')
    plt.ylabel('Valor', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3, axis='y')
    
    # Valores en las barras
    ax = plt.gca()
    for i, metodo in enumerate(metodos_estrategias):
        if metodo in trading_results:
            resultado = trading_results[metodo]
            valores = [
                resultado['total_return'],
                resultado['sharpe_ratio'],
                resultado['win_rate'] * 100
            ]
            for j, valor in enumerate(valores):
                ax.text(j + i * ancho_barra, valor + 0.01, f'{valor:.3f}', 
                       ha='center', va='bottom', fontsize=10)
    
    # Subgráfica 3: Betas comparadas
    plt.subplot(2, 2, 3)
    
    for metodo, resultado in trading_results.items():
        beta_vals = resultado['beta']
        series = range(1, 6)
        plt.plot(series, beta_vals[:5],
                color=colores[metodo], marker='o', linewidth=2,
                label=metodo, markersize=6)
    
    plt.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Beta1 Verdadero')
    plt.axhline(y=-1, color='gray', linestyle='--', alpha=0.5, label='Beta2 Verdadero')
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.title('COMPARACIÓN DE BETAS (Primeras 5 Series)', fontsize=14, fontweight='bold')
    plt.xlabel('Número de Serie', fontsize=12)
    plt.ylabel('Valor de Beta', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Subgráfica 4: Retornos por trade
    plt.subplot(2, 2, 4)
    
    for metodo, resultado in trading_results.items():
        retornos = resultado['trading_returns']
        if len(retornos) > 0:
            plt.plot(range(len(retornos)), retornos,
                    color=colores[metodo], marker='o', linewidth=1,
                    label=metodo, alpha=0.7, markersize=4)
    
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    plt.title('RETORNOS INDIVIDUALES POR TRADE', fontsize=14, fontweight='bold')
    plt.xlabel('Número de Trade', fontsize=12)
    plt.ylabel('Retorno del Trade', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('retornos_3_estrategias_T200.png', dpi=300, bbox_inches='tight')
    plt.show()
