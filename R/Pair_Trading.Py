# pair_trading.py
import numpy as np
import pandas as pd
from typing import Tuple, Dict, Optional
import warnings

def estrategia_pair_trading(
    X_train: np.ndarray,
    X_test: np.ndarray,
    beta: np.ndarray,
    config: Optional[Dict] = None
) -> Dict:
    """
    Implementa estrategia de pair trading basada en cointegración
    
    Parámetros:
    -----------
    X_train : np.ndarray
        Datos de entrenamiento para calcular estadísticas, forma (T_train, N)
    X_test : np.ndarray
        Datos de test para trading, forma (T_test, N)
    beta : np.ndarray
        Vector de cointegración, forma (N,) o (N, r)
    config : dict, opcional
        Configuración de la estrategia con keys:
        - entry_threshold: Z-score para entrada (default 1.5)
        - exit_threshold: Z-score para salida (default 0.5)
        - stop_loss: Z-score para stop loss (default None)
        - take_profit: Z-score para take profit (default None)
        - max_position: posición máxima (default 1)
        - transaction_cost: costo por transacción (default 0)
    
    Retorna:
    --------
    resultados : dict
        Diccionario con métricas y series de trading
    """
    # Configuración por defecto
    default_config = {
        'entry_threshold': 1.5,
        'exit_threshold': 0.5,
        'stop_loss': None,
        'take_profit': None,
        'max_position': 1.0,
        'transaction_cost': 0.0,
        'verbose': False
    }
    
    if config is not None:
        default_config.update(config)
    
    cfg = default_config
    
    try:
        # Asegurar que beta sea vector 1D
        if len(beta.shape) > 1:
            if beta.shape[1] > 1:
                warnings.warn("Beta tiene múltiples columnas, usando primera columna")
            beta = beta[:, 0] if beta.shape[1] > 0 else beta.flatten()
        
        # 1. Calcular spread en train y test
        spread_train = X_train @ beta
        spread_test = X_test @ beta
        
        # 2. Estadísticas del spread en train
        mean_spread = np.mean(spread_train)
        std_spread = np.std(spread_train)
        
        if std_spread == 0:
            warnings.warn("Desviación estándar del spread es 0")
            std_spread = 1e-10
        
        # 3. Calcular Z-scores en test
        z_scores = (spread_test - mean_spread) / std_spread
        
        # 4. Generar señales de trading
        n_test = len(spread_test)
        position = np.zeros(n_test)  # Posición actual
        entry_price = np.zeros(n_test)  # Precio de entrada
        trades = []  # Lista de trades
        
        for t in range(1, n_test):
            current_z = z_scores[t]
            prev_position = position[t-1]
            
            # REGLAS DE TRADING
            
            # Si no tenemos posición
            if prev_position == 0:
                # Señal LONG: Z-score muy bajo
                if current_z <= -cfg['entry_threshold']:
                    position[t] = cfg['max_position']
                    entry_price[t] = spread_test[t]
                    trades.append({
                        'timestamp': t,
                        'type': 'LONG',
                        'entry_price': spread_test[t],
                        'entry_z': current_z
                    })
                
                # Señal SHORT: Z-score muy alto
                elif current_z >= cfg['entry_threshold']:
                    position[t] = -cfg['max_position']
                    entry_price[t] = spread_test[t]
                    trades.append({
                        'timestamp': t,
                        'type': 'SHORT',
                        'entry_price': spread_test[t],
                        'entry_z': current_z
                    })
                
                else:
                    position[t] = 0
            
            # Si tenemos posición LONG
            elif prev_position > 0:
                # Salir si Z-score vuelve a la media
                if current_z >= -cfg['exit_threshold']:
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                
                # Stop loss (si está configurado)
                elif (cfg['stop_loss'] is not None and 
                      current_z <= -cfg['stop_loss']):
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                    trades[-1]['stop_loss'] = True
                
                # Take profit (si está configurado)
                elif (cfg['take_profit'] is not None and 
                      current_z >= cfg['take_profit']):
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                    trades[-1]['take_profit'] = True
                
                else:
                    position[t] = prev_position
            
            # Si tenemos posición SHORT
            elif prev_position < 0:
                # Salir si Z-score vuelve a la media
                if current_z <= cfg['exit_threshold']:
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                
                # Stop loss (si está configurado)
                elif (cfg['stop_loss'] is not None and 
                      current_z >= cfg['stop_loss']):
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                    trades[-1]['stop_loss'] = True
                
                # Take profit (si está configurado)
                elif (cfg['take_profit'] is not None and 
                      current_z <= cfg['take_profit']):
                    position[t] = 0
                    trades[-1]['exit_price'] = spread_test[t]
                    trades[-1]['exit_z'] = current_z
                    trades[-1]['duration'] = t - trades[-1]['timestamp']
                    trades[-1]['take_profit'] = True
                
                else:
                    position[t] = prev_position
            
            # Copiar entry_price si mantenemos posición
            if position[t] != 0 and entry_price[t] == 0:
                entry_price[t] = entry_price[t-1]
        
        # 5. Calcular retornos
        spread_returns = np.diff(spread_test)
        trading_returns = position[:-1] * spread_returns
        
        # Aplicar costos de transacción
        position_changes = np.diff(position) != 0
        transaction_costs = np.abs(position_changes) * cfg['transaction_cost']
        trading_returns = trading_returns - transaction_costs[:-1]
        
        # 6. Calcular métricas
        total_return = np.sum(trading_returns)
        cumulative_returns = np.cumsum(trading_returns)
        
        # Sharpe Ratio (anualizado si hay suficientes datos)
        if len(trading_returns) > 1 and np.std(trading_returns) > 0:
            sharpe_ratio = (np.mean(trading_returns) / np.std(trading_returns)) * np.sqrt(252)
        else:
            sharpe_ratio = 0
        
        # Win rate y otras métricas de trades
        winning_trades = 0
        total_trades_closed = 0
        
        for trade in trades:
            if 'exit_price' in trade:
                total_trades_closed += 1
                if trade['type'] == 'LONG':
                    profit = trade['exit_price'] - trade['entry_price']
                else:  # SHORT
                    profit = trade['entry_price'] - trade['exit_price']
                
                if profit > 0:
                    winning_trades += 1
        
        win_rate = winning_trades / total_trades_closed if total_trades_closed > 0 else 0
        
        # 7. Preparar resultados
        resultados = {
            # Datos básicos
            'spread_train': spread_train,
            'spread_test': spread_test,
            'z_scores': z_scores,
            'position': position,
            'entry_price': entry_price,
            
            # Retornos
            'trading_returns': trading_returns,
            'cumulative_returns': cumulative_returns,
            'total_return': total_return,
            
            # Métricas de riesgo
            'sharpe_ratio': sharpe_ratio,
            'volatility': np.std(trading_returns) if len(trading_returns) > 0 else 0,
            'max_drawdown': self._calcular_max_drawdown(cumulative_returns),
            
            # Métricas de trading
            'win_rate': win_rate,
            'num_trades': len(trades),
            'num_trades_closed': total_trades_closed,
            'avg_trade_duration': np.mean([t['duration'] for t in trades if 'duration' in t]) if total_trades_closed > 0 else 0,
            'trades': trades,
            
            # Configuración
            'config': cfg,
            'beta_used': beta,
            'train_stats': {
                'mean': mean_spread,
                'std': std_spread
            }
        }
        
        if cfg['verbose']:
            print(f"Estrategia completada:")
            print(f"  Retorno total: {total_return:.4f}")
            print(f"  Sharpe Ratio: {sharpe_ratio:.4f}")
            print(f"  Win Rate: {win_rate*100:.2f}%")
            print(f"  Número de trades: {len(trades)}")
        
        return resultados
        
    except Exception as e:
        warnings.warn(f"Error en estrategia de pair trading: {str(e)}")
        return {}


def _calcular_max_drawdown(self, cumulative_returns: np.ndarray) -> float:
    """Calcula el máximo drawdown de una serie de retornos acumulados"""
    if len(cumulative_returns) == 0:
        return 0
    
    # Calcular valor acumulado
    running_max = np.maximum.accumulate(cumulative_returns)
    drawdowns = (cumulative_returns - running_max) / (running_max + 1e-10)
    
    return np.min(drawdowns) if len(drawdowns) > 0 else 0


def backtest_pair_trading(
    X_full: np.ndarray,
    beta_func,
    train_size: int = 100,
    test_size: int = 50,
    window_type: str = 'rolling',
    **kwargs
) -> Dict:
    """
    Backtesting completo de estrategia de pair trading
    
    Parámetros:
    -----------
    X_full : np.ndarray
        Datos completos
    beta_func : callable
        Función que estima beta (debe aceptar X y retornar beta)
    train_size : int
        Tamaño de ventana de entrenamiento
    test_size : int
        Tamaño de ventana de test
    window_type : str
        'rolling' o 'expanding'
    **kwargs : dict
        Argumentos adicionales para beta_func y estrategia_pair_trading
    
    Retorna:
    --------
    backtest_results : dict
        Resultados del backtest
    """
    T_total = len(X_full)
    all_results = []
    
    if window_type == 'rolling':
        step_size = test_size
    else:  # expanding
        step_size = test_size
    
    for start in range(0, T_total - train_size - test_size + 1, step_size):
        train_end = start + train_size
        test_end = train_end + test_size
        
        X_train = X_full[start:train_end]
        X_test = X_full[train_end:test_end]
        
        # Estimar beta en train
        beta = beta_func(X_train, **kwargs)
        
        if beta is not None:
            # Ejecutar estrategia en test
            resultados = estrategia_pair_trading(
                X_train, X_test, beta, **kwargs
            )
            
            if resultados:
                resultados['period'] = {
                    'train_start': start,
                    'train_end': train_end,
                    'test_start': train_end,
                    'test_end': test_end
                }
                all_results.append(resultados)
    
    # Consolidar resultados
    if not all_results:
        return {}
    
    backtest_results = {
        'period_results': all_results,
        'total_return': np.sum([r['total_return'] for r in all_results]),
        'avg_sharpe': np.mean([r['sharpe_ratio'] for r in all_results]),
        'avg_win_rate': np.mean([r['win_rate'] for r in all_results]),
        'total_trades': np.sum([r['num_trades'] for r in all_results]),
        'num_periods': len(all_results)
    }
    
    return backtest_results
